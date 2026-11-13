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

struct SwapEdgeSink;   // rcpgenerator/diameter_permute.hpp

struct PackingConfig {
    std::vector<double> box;
    // Legacy packer behavior: empty walls defaults to Ndim zeros before first use.
    std::vector<std::int8_t> walls;
    std::uint32_t neighbor_max = 0;
    std::uint32_t seed = 0;
    bool fix_height = false;
    // Diameter-SWAP (additive). swap == 0 -> the packing loop is byte-identical to
    // the pristine build; every field below is inert unless swap > 0. Reachable from
    // Packing.pack() because it travels in PackingConfig, not PackingRunOptions.
    double swap = 0.0;             // strength: candidate pairs per particle per sweep
    int    swap_interval = 1;      // steps between sweeps
    int    swap_mode = 5;          // 5 directed, relative-D window (default) | 4 directed rank window | 1 rank window | 3 systematic | 0 random | 2 neighbour
    int    swap_rank_window = 8;   // mode 1 partner rank offset in [1, window]
    double swap_size_tol = 0.2;    // modes 0/2: |dD| <= tol * <D>
    double swap_eps = 0.0;         // accept if dE < -eps * (mean contact energy)
    long   swap_start = 0;         // first step sweeps may run
    long   swap_stop = -1;         // last step (-1: none)
    int    swap_phase = 0;         // bitmask over mu_flag: 1 rung1, 2 rung0, 4 rung-1; 0 all
    int    swap_adam = 0;          // 0 keep, 1 zero m, 2 zero m and v for swapped particles
    int    swap_inject_adam = -1;  // [206] Adam reset for swap_inject participants ONLY; -1 = follow swap_adam
    int    swap_inject_freeze = 0; // [207] local quench: for K steps after an injection only the participants move
                                   //       (positions AND Adam moments of everyone else are held; kappa held too)
    int    swap_inject_accept = 0; // [207] 1 = at the end of the quench keep the exchange only if U <= U0*(1+tol),
                                   //       else REVERT (diameters swapped back, participants' x/m/v restored)
    double swap_inject_tol = 0.0;  // [207] tolerance for the acceptance test
    int    swap_inject_freeze_shell = 1; // [207] free set = participants + this many neighbour shells (0 = participants only)
    int    swap_max_per_particle = 4;  // cap on accepted swaps per particle per sweep
    int    swap_refresh = 0;       // 1: force a neighbour refresh of swapped particles (not needed for small |dD|)
    int    swap_backoff = 64;      // max multiplier of swap_interval after empty sweeps (1 = fixed interval)
    int    swap_exact = 0;         // 0: O(1) proxy from per-particle (g, z) (default); 1: exact dE over near-pair lists
    double swap_window_rel = 0.016; // mode 5: partner within |D_j/D_i - 1| <= this (scale-free window)
    double swap_relayout_ratio = 1.3; // if a swapped particle's D0 grew past this ratio since the last
                                      // pair-layout, recompute layout (spatial reorder) next step
    int    swap_rank_space = 1;    // 1: rank-space (streaming) proxy sweep; 0: random particle draws
    double swap_idle_frac = 1e-2;  // a sweep accepting <= idle_frac*N exchanges counts as idle (back-off)
    // (a, R) outcome controller (opt-in; overrides fixed strength/window semantics)
    int    swap_ctrl = 0;          // 1: quota-based draws + adaptive window
    double swap_a = 0.3;           // target accepted exchanges per particle per sweep
    double swap_R = 0.20;          // target acceptance fraction
    double swap_eta = 0.2;         // window feedback gain
    double swap_wmin = 5e-4;       // window floor
    double swap_wmax = 0.10;       // window ceiling (initial window)
    double swap_smax = 32.0;       // draw cap per particle per sweep
    double swap_floor_budget = 1.0; // ctrl at window floor: draw budget x (a/R) pp; 0 = legacy cap-burn
    double swap_rescue = 0.0;      // rescue channel: exact-dE candidates pp per rescue pass (0 = off)
    double swap_rescue_w = 0.05;   // rescue partner window |Dj/Di - 1| <= w
    int    swap_rescue_interval = 4; // steps between rescue passes
    double swap_cycle3 = 0.0;      // 3-cycle trio candidates pp per rescue pass (0 = off)
    // ---- session survivability [log 043] -------------------------------------------
    std::string checkpoint_path;   // non-empty: periodic bitwise state checkpoints to
                                   // this file (atomic tmp+rename), AND resume from it
                                   // at pack() start if it exists (header-validated)
    double checkpoint_interval_s = 1800.0; // wall seconds between checkpoint writes
    long   max_steps = 0;          // >0: hard step cap (fail-soft; result.capped set).
                                   // 0 = schedule default.
    // ---- mechanics instrumentation (ME0 [log 054]; all default off) ---------------
    long   swap_event_log_stride = 0;   // per-trade event log sampling (0 = off)
    long   swap_event_log_cap = 1000000; // max logged events
    int    swap_event_log_geom = 0;      // ME0-iv: +cage/fact fields (row width 12)
    long   checkpoint_every_steps = 0;   // >0: archival write every k steps (path.STEP)
    // ---- mu-hold actuator (MECHANICS stage 2b [log 071]; labeled variant, default off) ----
    long   mu_hold_start = -1;           // first step of the hold (-1 = off)
    long   mu_hold_steps = 0;            // hold length in steps (0 = off)
    double mu_hold_mult = 1.0;           // mu_eff = mult * mu_1 (initial rung-1 mu) during the hold
    int    mu_hold_shift_deadline = 1;   // 1: a hold that starts in rung 1 delays the deadline by its length
    int    record_visibility = 0;        // 1: per-step visibility P(z>=1) and mu_eff histories
    // ---- visibility controller prototype (stage 2b Part 5 [log 073]; labeled, default off) ----
    // Rides on the mu_hold window: within [mu_hold_start, +mu_hold_steps) the
    // multiplier is adapted every mu_ctrl_interval steps to hold v = P(z>=1) at
    // mu_ctrl_v_target (deadband), clamped to [mu_ctrl_min_mult, mu_ctrl_max_mult].
    int    swap_mu_ctrl = 0;
    double mu_ctrl_v_target = 0.30;
    double mu_ctrl_deadband = 0.05;
    double mu_ctrl_eta = 0.5;            // mult <- mult * exp(eta * (target - v)) outside the deadband
    long   mu_ctrl_interval = 50;
    double mu_ctrl_min_mult = 0.1;
    double mu_ctrl_max_mult = 1.0;
    // ---- stage 3 [log 082]: geometric sort step (labeled, default off) ----
    int    swap_geo = 0;                 // 1 replace sweep on geo steps; 2 hybrid (both)
    double swap_geo_T = 0.0;
    double swap_geo_T_end = -1.0;        // [135] H9 anneal floor T (<0 = fixed swap_geo_T)
    double swap_geo_T_decay = 1.0;       // [135] H9 per-pass T multiplier (<1 anneals)
    double swap_geo_w = 0.05;
    long   swap_geo_B = 2000;
    long   swap_geo_k = 200;
    long   swap_geo_start = 0;
    long   swap_geo_stop = -1;
    long   swap_geo_reopen = 0;    // [101] max convergence reopens: at OSC DONE, force a
                                   // geo pass and re-enter FINAL_DESCENT (same constants,
                                   // fresh budget); stops early once a reopen trades nothing
    int    swap_geo_cycle = 0;     // [117] H2: apply 3-cycles after each geo pass (k-opt escape)
    int    swap_geo_cycle_cap = 12;// per-rank partner scan cap for the cycle search
    int    swap_geo_cycle4 = 0;    // [194] H34: 4-cycle pass after the 3-cycle pass (k-opt rung 4)
    int    swap_geo_cycle4_cap = 6;// per-rank partner scan cap for the quadruple search
    int    swap_geo_capk = 0;      // [121/122] H7: cage-capacity sort; K smallest per-nbr caps (0=off)
    int    swap_geo_capmode = 3;   // 1=K-th smallest, 2=mean-of-K, 3=tiered c_min-else-cage
    int    swap_geo_bottleneck = 0;// [146/148/152] H12: bottleneck-relief. 0=off; 1=terminal (reopen-loop final squeeze); 2=H12b under-overlap at the mu_1 rung (mu_flag==1, additive)
    int    swap_geo_cage = 0;      // [129/130/131] H8: true-cage overlap-energy pairwise move (0=off,1=on)
    int    swap_geo_cage_move = 0; // [133] DIAGNOSTIC: also move swapped particle to its cage centre (0=off; position move, OUTSIDE C5 - insight only)
    int    swap_geo_nonsep = 0;    // [172] H19: non-separable comparator rescore of geo candidates (0=off,1=on; C5-clean, diameters-only)
    double swap_geo_ctarget = 0.5; // [172] H19: concave clearance cap in units of mean radius
    int    swap_geo_nsk = 4;       // [172] H19: top-K multiplier over geo_B for the non-separable rescore set
    long   swap_geo_nonsep_start = 0;   // [190] H26: coupled rescore (nonsep==1) fires only for step >= start
    long   swap_geo_nonsep_stop = -1;   // [190] H26: ... and step < stop (-1 = open)
    long   swap_geo_cage_start = 0;// [131] H8 HANDOFF: step at which geo yields to the cage pass (0=both)
    long   swap_geo_capk_start = 0;// [127] H7 HANDOFF: step at which geo yields to cage. Before this
                                   // step geo_pass runs (c_min sort during compression); from this
                                   // step geo_pass is OFF and capk_pass runs ALONE on the geo cadence
                                   // (cage polish of the near-jammed residual). 0 = capk from geo_start
                                   // (the old competing-pass behavior). Avoids the geo<->cage thrash [126].
    int    swap_sweeps_off = 0;    // 1: suppress normal sweeps (surgical mode)
    std::vector<double> swap_inject;    // flattened (step, i, j) injection triples
    // ---- [219] DISTRIBUTION MORPH (labeled, default off): evolve the base-diameter multiset D0 along the
    // rank-preserving quantile interpolation (1-lambda)*Q_start + lambda*Q_target, renormalised each update
    // so sum(D0^3) is constant; a feedback controller advances lambda while max overlap stays in band. ----
    std::vector<double> morph_target;   // N target relative diameters (any order; sorted internally); empty = off
    long   morph_start = 0;             // first step at which the morph may update (absolute step)
    long   morph_interval = 10;         // steps between updates (1 = continuous, rank order cached by epoch; 10 measured better at equal rate [220])
    double morph_kappa_damp = 1.0;      // multiplier on kappa's Adam step while the morph is advancing (lambda < 1):
                                        //   < 1 lets positions, not the global scale, absorb the transient imbalance
    double morph_band_floor = 5e-4;     // the band reference is max(max overlap at arming, this floor) - absolute overlap units
    double morph_margin = 0.5;          // the increment may grow only while max overlap < margin * band * reference
    long   morph_max_hold_streak = 2000;// consecutive held updates that declare a kappa-lag STALL: morphing stops, run ends
    double morph_kappa_drop = 1e-3;     // deflation guard: hold while kappa < baseline * (1 - this) - the kappa-lag signal itself
    long   morph_kappa_tau = 100;       // updates over which the baseline tracks kappa (a quasi-static drift passes; a burst trips)
    double morph_rate0 = 1e-3;          // initial lambda increment per update
    double morph_band = 2.0;            // hold when max overlap > band * reference (measured at morph_start)
    double morph_grow = 1.25;           // increment growth factor after morph_grow_after in-band updates
    double morph_shrink = 0.7;          // increment shrink factor on a hold (0.5 ratcheted stiff states to the floor [220])
    double morph_rate_min = 1e-5;       // floor of the increment: below it a hold only waits (costs steps, not progress)
    long   morph_grow_after = 5;        // consecutive in-band updates before growing the increment
    double morph_rate_max = 3e-5;       // cap on the lambda increment: the rate at which the S 4->3 read-out is rate-independent
                                        //   ([220] plateau); feedback alone settles at 3-5e-4 and reads 0.005 low
    double morph_settle_flat = 0.0;     // >0: after lambda reaches 1, END the run once max overlap is in band and the
                                        //     phi range over morph_settle_window steps is below this value (settled read-out)
    long   morph_settle_window = 1000;  // window (steps) for the settle test
    // ---- [219] mu pin (labeled, default off): mu fed to the forces is held at this absolute value and the
    // schedule's rung/decay/oscillation transitions are suspended (the run continues to max_steps). ----
    double mu_pin = 0.0;                // <= 0 = off
    // ---- [221] QUENCH to force balance (labeled, default off): FIRE position updates replace Adam, kappa frozen,
    // schedule suspended; ends at the per-particle residual tolerance on the load-bearing backbone (z >= Ndim + 1,
    // wall contacts counted) or at a give-up trigger. ----
    int    quench_mode = 0;             // 0 off, 1 FIRE
    double quench_tol = 1e-10;          // max over backbone of |F_i| / sum_j |f_ij|
    long   quench_max_iter = 200000;    // iteration cap (steps)
    double quench_dt0 = 0.05, quench_dt_max = 0.5, quench_alpha0 = 0.1;
    int    quench_mass_scaled = 0;      // 1: acceleration = F / m_i with m_i = z_i / median(z) (preconditioned FIRE)
    long   quench_report_every = 500;   // quench_log stride
    int quench_two_level = 0;           // [223] 1: FIRE over the heavy set only; light particles (below quench_light_frac of the typical load,
                                        //       rattlers included) relax by an overdamped local step x += omega F/(K z) against the frozen heavy set
    double quench_light_frac = 0.3;     // load fraction (of 4 x mean contact force) below which a particle is "light"
    double quench_light_omega = 1.0;    // relaxation factor of the light step
    double quench_freeze_damp = 0.0;    // [224] activity-gated dissipation: velocity damping factor per step on quiet particles (0 = off)
    double quench_freeze_frac = 0.1;    // [224] quiet = EMA speed below this fraction of the median
    long quench_warm_steps = 0;         // [224] controlled warm phase: iterations with velocities rescaled and no P<0 resets (0 = off)
    double quench_warm_amp = 0.1;       // [224] warm-phase mean displacement per step, in mean overlaps
    int quench_class_every = 100;       // iterations between heavy/light reclassifications
    int quench_mix_mode = 0;            // FIRE mixing: 0 global norms (Bitzek), 1 per-particle norms, 2 per-particle norms + per-particle P<0 reset [222]
    double quench_kappa_mult = 1.0;     // kappa <- kappa * mult at the first quench iteration (quench at a chosen pressure)
    std::vector<double> quench_held;    // indices of held particles (masked out; their net force is reported)
    // [225] give-up triggers (each 0 = off; the result reports which one ended the quench in quench_status)
    long   quench_stall_window = 0;     // stop as "stalled" when the best residual has not fallen by quench_stall_factor
    double quench_stall_factor = 0.5;   //   within this many iterations
    double quench_dt_min = 0.0;         // stop as "step collapse" when FIRE's time step falls below this
    // Ellipsoids (additive; rcpgenerator/ellipsoid.hpp). Empty -> spheres, the code path is untouched.
    // Ndim-1 numbers in decreasing order, 1 >= ratio_2 >= ... > 0: axis k (k = 2..Ndim) has length D x ratio_k;
    // D is the first and LONGEST axis. One shape for every particle; for a shape per particle use
    // PackingInput::aspect_ratios instead (not both). Must come with PackingInput::orientations.
    std::vector<double> aspect_ratios;
    int ell_iters = 5;                  // fixed number of contact-function iterations per overlapping pair
    double ell_rot_scale = 1.0;         // orientation step = ell_rot_scale x alpha / (bounding radius)
    double ell_track_tol = 1e-6;        // lam tracking tolerance inside a run (0: fixed ell_iters + 1 evaluations)
    // Diameter permutation with ellipsoids (BUILD.md section 8; see ell_permute_pass in ellipsoid.hpp)
    // Rods (additive; rcpgenerator/rod.hpp). rod_capped applies when PackingInput::rod_length / rod_radius are set:
    // 1 capped (half-sphere ends), 0 uncapped (flat ends).
    int rod_capped = 1;
    // Rounded rods (rod.hpp, "rounded rods"): relative rounding eps of the rod's edges. rod_round applies to single
    // evaluations (evaluate_once); quench_rod_round to a quench run. 0: the sharp rod of the packing loop.
    double rod_round = 0.0;
    double quench_rod_round = 0.0;
    // Trajectory samples spaced evenly in log(step) instead of every trajectory_interval steps: this many per
    // decade of steps (0: off, the interval is used). Only read when positions are captured.
    int trajectory_per_decade = 0;
    int ell_permute_orient = 0;         // 0 (default, Ken 2026-11-12 on the D5/D6 table): diameters and shapes are
                                        //    exchanged, orientations stay with the sites; 1: orientations travel
                                        //    with the particles (the move of BUILD.md section 8, kept as an option);
                                        // 2: diagnostic: decided as in 0, moved as in 1
    int ell_permute_partners = 64;      // partners examined per particle per pass (the sphere engine's cap)
    int ell_permute_guard = 1;          // 1: an exchange must also not raise the two particles' cage overlap
                                        //    (not applied to a pair of spheres); 0: room gain alone
    long ell_capture = 0;               // test hook: n > 0 records the state and the engine's forces, torques and
                                        // energy at the force pass that follows the n-th pass with an exchange
};

struct PackingInput {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
    // Cascade (additive): per-particle flag, 1 = diameter frozen (exempt from the
    // global kappa scaling). Empty / all-zero -> standard one-shot path (unchanged).
    std::vector<std::int8_t> fixed;
    // Ellipsoids (additive): one rotation matrix per particle, row-major Ndim x Ndim (Ndim*Ndim numbers),
    // columns = body axes in the lab frame. Empty -> spheres. Both-or-neither with the aspect ratios.
    std::vector<std::vector<double>> orientations;
    // Ellipsoids, a shape per particle (additive): N rows of Ndim-1 ratios, each obeying the decreasing-order
    // rule. Empty -> the single shape of PackingConfig::aspect_ratios (or spheres); supplying both is an error.
    // A row follows its particle wherever the diameter and the orientation go.
    std::vector<std::vector<double>> aspect_ratios;
    // Rods (additive): straight length and radius in units of the particle's D, each one number (every particle) or
    // N numbers. Both or neither; with orientations; never together with aspect ratios.
    std::vector<double> rod_length, rod_radius;
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
    std::vector<float> zpos_history;          // record_visibility: fraction of particles with z >= 1
    std::vector<float> mu_eff_history;        // record_visibility: mu fed to forces (hold-aware)
    long mu_hold_applied = 0;                 // steps on which the hold overrode mu
    double mu_ctrl_final_mult = 0.0;          // controller multiplier at the end of the window
    long swap_geo_trades = 0;                 // accepted geometric-pass exchanges
    long swap_geo_passes = 0;
    long swap_geo_reopens = 0;                // [101] convergence reopens actually taken
    long swap_geo_cycle_trades = 0;           // [117] accepted 3-cycle rotations
    long swap_geo_cycle4_trades = 0;          // [194] accepted 4-cycle rotations
    long swap_geo_capk_trades = 0;            // [121/122] accepted cage-capacity swaps
    long swap_geo_cage_trades = 0;            // [129/130/131] accepted overlap-energy pairwise swaps
    long swap_geo_bottleneck_trades = 0;      // [146/148] accepted bottleneck-relief swaps
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
    // [219] distribution morph diagnostics: 6 doubles per update [step, lambda, increment, max_overlap, held, sum D0^3]
    std::vector<double> morph_log;
    double morph_lambda = 0.0;                // final morph fraction (1.0 = target reached)
    long   morph_updates = 0;                 // updates applied (lambda advanced)
    long   morph_holds = 0;                   // updates held (out of band)
    long   morph_stalled = 0;                 // 1 if the run ended on a kappa-lag stall (hold streak > morph_max_hold_streak)
    double morph_rate_ceiling = 0.0;          // final ceiling of the increment (lowered on every band trip)
    // [221] quench diagnostics
    double quench_residual = 0.0;             // final max backbone residual
    double quench_global_residual = 0.0;      // final RMS net force / mean contact force
    long   quench_iters = 0;
    int    quench_converged = 0;
    int    quench_status = 0;                 // [225] 0 not run, 1 converged, 2 max_iter, 3 stalled, 4 unjammed,
                                              //       5 step collapse, 6 non-finite force or energy
    std::vector<double> quench_log;           // 6 per report: iter, residual, global, U, n_backbone, dt
    std::vector<double> quench_held_force;    // Ndim doubles: net contact force on the held set at exit
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
    // Diameter-SWAP: number of extra spatial reorders triggered by the capacity
    // guard (0 when swap == 0).
    std::size_t swap_relayouts = 0;
    // rescue channel totals (proxy-vs-exact cross-table; [log 034])
    std::size_t swap_rescue_candidates = 0;
    std::size_t swap_rescue_accepted = 0;
    std::size_t swap_rescue_false_rej = 0;
    std::size_t swap_rescue_false_acc = 0;
    std::size_t swap_cycle3_candidates = 0;
    std::size_t swap_cycle3_accepted = 0;
    std::size_t swap_cycle3_blocked = 0;   // accepted cycles with all pair swaps uphill
    bool capped = false;                   // run ended at the max_steps cap (fail-soft)
    std::size_t resumed_from_step = 0;     // 0 = fresh run; else checkpoint resume point
    std::size_t checkpoint_last_step = 0;  // step of the last successful checkpoint write
    std::vector<double> swap_events;
    int swap_event_width = 8;          // doubles per event row       // per-trade event log (8 doubles/event)
    // Ellipsoids (additive): final orientations, in the same particle order as positions; empty for spheres.
    std::vector<std::vector<double>> orientations;
    // Per-particle shapes in the same particle order as positions; filled only when the input carried a table.
    std::vector<std::vector<double>> aspect_ratios;
    std::vector<double> rod_length, rod_radius;   // rods: per particle, in the output particle order
    // Permutation with ellipsoids: totals, and the ell_capture record (empty unless requested and reached).
    std::size_t ell_permute_candidates = 0;
    long ell_capture_step = -1;
    double ell_capture_energy = 0.0;
    std::vector<double> ell_capture_positions, ell_capture_diameters, ell_capture_orientations,
                        ell_capture_ratios, ell_capture_forces, ell_capture_torques, ell_capture_exchanges,
                        ell_capture_z;
    double torque_magnitude = 0.0;     // mean |generalised torque| / (mean contact force x mean radius) at the end
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
    std::vector<std::vector<std::vector<double>>> orientations;   // ellipsoids only; one entry per sample
    // Which particle sits in which slot of each sample: origins[f][k] = the index at the start of the run of the
    // particle in slot k of sample f (the loop renumbers particles for cache locality while it runs). One entry
    // per sample that has positions.
    std::vector<std::vector<std::size_t>> origins;
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

// ---- single evaluations (additive; used by the ellipsoid tests and for reporting) --------------------------
// One evaluation of forces, torques, energy and the kappa gradient on the configuration as given (no rescale,
// no steps), through the same input rule as a run: no aspect ratios and no orientations -> the sphere kernel
// get_forces_nd_3; both -> the ellipsoid kernel; one without the other -> error. All pairs are listed.
struct SingleEvaluation {
    std::vector<std::vector<double>> forces;    // N x Ndim, force on each particle (-dE/dx)
    std::vector<std::vector<double>> torques;   // N x Ndim(Ndim-1)/2 (zeros for spheres)
    double U = 0.0;             // the code's diagnostic (square of the mean fractional overlap)
    double dkappa = 0.0;        // the code's kappa gradient, including the mu term
    double max_overlap = 0.0;   // max over pairs of 1 - s  (1 - dist/r_ij for spheres)
    double energy = 0.0;        // sum over contacts of (K/2)(sigma - d)^2 (+ wall terms)
    double Fmean = 0.0;
    double Lc = 0.0;
    std::vector<std::size_t> z;
};
SingleEvaluation evaluate_once(const struct PackingInput& input, const struct PackingConfig& config, double mu);

// Exact shrink-to-contact of a configuration: the largest uniform factor (<= 1) by which every particle can be
// scaled about its centre so that no pair (and no flat wall) overlaps, computed to convergence. Spheres too.
struct ContactScale {
    double scale = 1.0;             // min(1, min over contacts of s)
    std::size_t pair_i = 0, pair_j = 0;   // the tightest pair (pair_j == N: a wall)
    std::size_t n_overlapping = 0;  // contacts with s < 1 before the shrink
    double particle_volume = 0.0;   // total particle volume before the shrink
};
ContactScale contact_scale(const struct PackingInput& input, const struct PackingConfig& config);

// Ellipsoid permutation outside a run (tests and inspection). The state is the input's; all pairs are listed.
struct PermuteOnce {
    std::vector<double> diameters;
    std::vector<std::vector<double>> aspect_ratios;   // N rows
    std::vector<std::vector<double>> orientations;    // N flat row-major matrices
    std::vector<double> rod_length, rod_radius;       // rods: per particle after the pass
    std::vector<double> exchanges;                    // (i, j, gain) per accepted exchange
    std::size_t candidates = 0;
};
// One permutation pass on the configuration as given: window config.swap_geo_w, at most config.swap_geo_B
// exchanges (<= 0: no limit), config.ell_permute_partners (<= 0: every partner in the window),
// config.ell_permute_orient.
PermuteOnce permute_once(const struct PackingInput& input, const struct PackingConfig& config);
// One forced exchange of particles i and j (the move itself, no decision).
PermuteOnce exchange_once(const struct PackingInput& input, const struct PackingConfig& config, std::size_t i,
                          std::size_t j);

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
    double kappa = 1.0,
    // Cascade (additive): per-particle kappa-radius (fixed->0) + active-only diameter
    // sum. Defaults nullptr/-1 -> standard path, byte-identical (masked dkappa terms
    // reduce exactly to the originals).
    const double* rk_data = nullptr,
    double sumD_active = -1.0,
    // Diameter-SWAP (additive): when non-null, (a) if sink->want_edges, every
    // candidate pair with d2 < ((Di+Dj)/2 + sink->margin)^2 is appended to
    // sink->per_thread[tid] as (i, j, d2); (b) if sink->want_gsum, the absolute
    // overlap of every accepted contact is accumulated per particle into
    // sink->gsum. nullptr (the default, and always when swap == 0) -> nothing,
    // byte-identical.
    struct SwapEdgeSink* swap_sink = nullptr);

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
