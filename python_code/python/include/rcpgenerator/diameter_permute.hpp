#pragma once
// Diameter permutation: permutes which particle carries which diameter. Additive to the
// packing loop; nothing here runs unless it is switched on (internal fields keep a legacy swap_ prefix).
//
// A "sweep" runs at fixed positions, right after the force evaluation. Default
// (exact = 0): each candidate pair (i, j) is judged by the O(1) first-order estimate
//     dE ~ (g_grow - g_shrink) dR + 1/2 (z_grow + z_shrink) dR^2
// of the change of the solver's objective E = 1/2 sum max(0, (D_i+D_j)/2 - r_ij)^2,
// where g is a particle's summed absolute overlap and z its contact count, both
// accumulated by the force loop. Candidates are the particle's rank neighbours in
// diameter (mode 5: within a relative window; the direction -- grow or shrink -- is
// chosen by the particle's load), organised in rank space so the sweep streams.
// Research path (exact = 1): exact dE over near-pair edges collected by the force
// loop. Only D0/D are mutated; the caller handles cache invalidation.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace rcpgenerator {

struct SwapParams {
    double swap = 0.0;          // strength: candidate pairs per particle per sweep. 0 = off.
    int    interval = 1;        // steps between sweeps
    int    mode = 5;            // 5 directed, relative-D window (default) | 4 directed rank window | 1 rank window | 3 systematic | 0 random | 2 neighbour
    int    rank_window = 8;     // mode 1: partner rank offset uniform in [1, rank_window]
    double size_tol = 0.2;      // modes 0/2: reject |D_i - D_j| > size_tol * <D>
    double eps = 0.0;           // accept if dE < -eps * e_scale (e_scale = mean contact energy)
    long   start = 0;           // first step at which sweeps may run
    long   stop = -1;           // last step (-1 = no limit)
    int    phase = 0;           // bitmask over mu_flag: 1 rung1, 2 rung0, 4 rung-1; 0 = all
    int    adam = 0;            // 0 keep optimizer state; 1 zero m; 2 zero m and v (swapped particles)
    int    max_per_particle = 4; // accepted swaps per particle per sweep
    double window_rel = 0.016;  // mode 5: partner within |D_j/D_i - 1| <= window_rel (scale-free)
    int    exact = 0;           // 0: O(1) proxy from (g, z) (default); 1: exact dE over near-pair lists
    int    rank_space = 1;      // proxy sweep organised in rank space (streaming; default) vs particle draws
    // ---- (a, R) outcome controller (opt-in; requires mode 5, rank_space, proxy) ----
    int    ctrl = 0;            // 1: quota-based draws + adaptive window
    double target_a = 0.3;      // target ACCEPTED exchanges per particle per sweep
    double target_R = 0.20;     // target acceptance fraction (accepted/candidates)
    double ctrl_eta = 0.2;      // w <- w * exp(eta * (A - R)), clamped
    double w_min = 5e-4;        // window floor (starvation at floor = lock-in signal)
    double w_max = 0.10;        // window ceiling (also the initial w)
    double s_max = 32.0;        // draw cap per particle per sweep (cost bound)
    double floor_budget = 1.0;  // ctrl, window at w_min: draw budget factor x (target_a/target_R)
                                // per particle (healthy-phase budget) instead of draw-to-cap;
                                // 0 = legacy cap-burn. s_max stays an upper bound.
    double idle_frac = 1e-2;    // starvation threshold (shared with the fixed-w mode)
    // ---- rescue channel (opt-in; requires mode 5; [log 034]) ----------------------
    // A rare, serial, EXACT-dE candidate stream at a wide relative-D window, run
    // after the main proxy sweep on rescue steps. Purpose: large-delta exchanges the
    // O(1) proxy would falsely reject (the proxy is exact only for a fixed contact
    // set). Instrumented: proxy-vs-exact decision disagreement is counted.
    double rescue = 0.0;        // candidates per particle per rescue pass. 0 = off.
    double rescue_w = 0.05;     // partner window |D_j/D_i - 1| <= rescue_w
    int    rescue_interval = 4; // steps between rescue passes (edge collection cost
                                // is paid only on these steps)
    // ---- 3-cycles (opt-in; same pass/window/interval as the rescue channel) -------
    // Diameter rotations i<-D_j, j<-D_k, k<-D_i (both orientations tried), all three
    // mutually within rescue_w, decided by exact dE. Escapes states where every PAIR
    // exchange is uphill; accepted cycles whose three pairwise swaps are all uphill
    // are counted (the pairwise-blocked pool, [log 035]).
    double cycle3 = 0.0;        // trio candidates per particle per rescue pass. 0 = off.
    // ---- per-trade event log (ME0-i [log 054]; default off) -----------------------
    // Deterministic sampling (per-band counter, NO RNG consumption): every
    // stride-th accepted trade appends 8 doubles to event_buf_:
    //   [step, i_grow, j_shrink, delta_rel, delta_abs(dR), relief(-dE; proxy value in the
    //    banded channel, exact in rescue/cycle; NaN for injections),
    //    z_grow, z_shrink]  (z = pre-trade contact counts; -1 if unknown)
    long   event_log_stride = 0;   // 0 = off; 1 = every accept; k = every k-th per band
    long   event_log_cap = 1000000; // max events kept (buffer stops, flag set)
    // ME0-iv [log 067]: geometry fields per event (row widens 8 -> 12):
    //   [.., cage_grow, cage_shrink, fact_grow, fact_shrink]
    // cage_i = mean over near-pair neighbours k of (r_ik - R_k)  (a radius);
    // fact_i = z_i / n_nbr_i (share of near-pair neighbours in overlap).
    // Needs near-pair edges on every sweep step (rcp_generator requests them
    // when this flag is on); NaN when unavailable. Default off.
    int    event_log_geom = 0;
    // ---- stage 3 [log 082]: geometric sort step on the inscribed radius ---------
    // geo=1: on geo steps the proxy sweep is REPLACED by a geo pass; geo=2: both.
    // A pass proposes in-window exchanges (|D_j/D_i - 1| <= geo_w), scores them by
    // the misfit change dM = 2 dR (c_j - c_i) with c = inscribed radius from the
    // near-pair edges, accepts with prob min(1, exp(-dM/T)) (T=0: greedy, sorted),
    // each particle once per pass, at most geo_B per pass, every geo_k steps.
    int    geo = 0;
    double geo_T = 0.0;
    double geo_w = 0.05;
    long   geo_B = 2000;
    long   geo_k = 200;
    long   geo_start = 0;
    long   geo_stop = -1;
    int    geo_max_partners = 64;   // per-rank partner scan cap (evenly strided)
    // ---- H9 [log 135]: SIMULATED ANNEALING of the geo swap -----------------------
    // geo_T is the acceptance temperature (0 = greedy 2-opt). To escape the greedy
    // local optimum, ANNEAL: start at geo_T and multiply by geo_T_decay each pass,
    // floored at geo_T_end, so early passes accept some worsening swaps (escape) and
    // late passes are exact (T->0). geo_T_decay=1 (or geo_T_end<0) => fixed T =
    // backward-compatible (parity). Only the swap SEARCH is affected, not the solver.
    double geo_T_end = -1.0;        // floor temperature; <0 disables annealing (fixed geo_T)
    double geo_T_decay = 1.0;       // per-pass multiplier on the live temperature (<1 anneals)
    // ---- H2 [log 107/117]: 3-cycle sort move, gated, runs AFTER geo_pass -------
    // geo_cycle=1: after the pairwise geo pass, apply beneficial 3-cycles that NO
    // transposition among the triple can reach (the k-opt escape of a 2-opt-stable
    // state). Objective J = sum D_i*cmin_i; a rotation's dJ is NOT a sum of pairwise
    // dJ. Only fires where all three pairwise moves are non-improving. Inert at 0.
    int    geo_cycle = 0;
    int    geo_cycle_cap = 12;      // per-rank partner scan cap for the cycle search
    // [194] H34: geo_cycle4=1 adds a 4-cycle pass AFTER the 3-cycle pass: rotations of
    // diameters among an in-window quadruple accepted only when improving AND better than
    // every lower-order move within it (6 transpositions, 3 double-transpositions, 8 3-cycles).
    int    geo_cycle4 = 0;
    int    geo_cycle4_cap = 6;      // per-rank partner scan cap for the quadruple search (cost ~N*cap^3)
    // ---- H7 [log 121/122]: cage-capacity sort move, gated, runs AFTER geo_pass ---
    // c_min = min over neighbours of (r_ij - R_j) is a single-facet estimator and is
    // degenerate (-> R) at a jammed terminal. 3D translational jamming needs >= 4
    // contacts, so the K SMALLEST per-neighbour capacities describe the cage. This
    // pass keeps cap[i][0..K-1] (K = geo_capk) and re-ranks D against a cage signal.
    // It is OVERLAP-TOLERANT by design (no feasibility reject) - the transient
    // overlap is resolved by the solver's relaxation into the cage room. Inert at 0.
    //   geo_capk    : 0 = off; K = how many smallest per-neighbour caps to keep (>=1).
    //   geo_capmode : 1 = K-th smallest cap[K-1]; 2 = mean of the K smallest;
    //                 3 = tiered "A else B": accept if the c_min(=cap[0]) match
    //                     improves, ELSE if the cage-composite match improves
    //                     (this taps the residual c_min-sorting is blind to).
    int    geo_capk = 0;
    int    geo_capmode = 3;
    // ---- H8 [log 129/130/131]: true-cage OVERLAP-ENERGY pairwise move, gated -------
    // The swap runs ABOVE jamming under compression, where the overlap IS the
    // relaxation signal (user). The overlap of a trial diameter D at site i with a
    // neighbour k is max(0, D/2 - gap_ik) with gap_ik = r_ik - R_k (independent of
    // i's own diameter, so frozen at pass start like geo's c_min). The cage overlap
    // ENERGY E(D@i) = sum_k max(0, D/2 - gap_ik)^2 over i's K nearest contacts. A
    // pair swap is accepted iff it lowers total energy: E(Dj@i)+E(Di@j) < keep.
    // This is the pairwise-COMPLETENESS certificate c_min cannot give: when
    // geo_cage_trades -> 0, no overlap-improving pairwise swap remains (2-opt on the
    // true objective). Static probe [130/130b]: captures ~all the supply c_min misses
    // (recentering to the tetra centre adds <2.5%, so scored at the CURRENT position;
    // no position move - user). OVERLAP-TOLERANT. Inert at 0.
    //   geo_cage : 0 = off; 1 = on (overlap-energy pairwise pass after geo_pass).
    int    geo_cage = 0;
    // H12 [146/148]: bottleneck-relief "final squeeze" pass, run in the reopen loop
    // (relieve tightest contacts via accommodating swaps; gentle small-mu relax
    // between passes = the reopen re-descend). 0 = off; 1 = on. Feasibility-respecting.
    int    geo_bottleneck = 0;
    // geo_cage_move : DIAGNOSTIC ONLY [133]. 0 = diameters-only (pure permutation,
    // valid under C5). 1 = ALSO move each swapped particle to its cage centre at
    // accept time (deepest point of its cell). This is a POSITION move - OUTSIDE the
    // C5 move space - so any phi it produces is an INSIGHT number (does the cage
    // "close higher" if the particle recentres?), never a milestone phi. Static
    // probe [133]: recentring attacks ~27% of overlap energy vs 2.3% for the swap,
    // so this measures whether that headroom is dynamically realizable. Inert at 0.
    int    geo_cage_move = 0;
    // ---- H19 [log 172]: NON-SEPARABLE comparator (the frontier lever) --------------
    // geo's objective sum_i D_i*cmin_i is SEPARABLE (cmin frozen; the swap gain counts
    // only i,j's own terms), so pairwise 2-opt == global (H11) and cyc/cage/bneck add ~0.
    // H19 rescoring: maximise J' = sum_i g(cmin_i) with cmin_i = min_k (r_ik - 0.5 D_k)
    // LIVE and g concave-capped g(c) = min(c, ctarget) (relieve the tightest cages first,
    // saturating). Swapping D_i<->D_j changes cmin of EVERY neighbour of i and j (the
    // third-party term the geo gain drops) => dJ' is non-separable. Two-stage: rank by the
    // cheap separable gain, rescore the top (nsk*geo_B) candidates by the true dJ' over
    // N(i) U N(j), accept greedily by dJ'. Diameters-only (C5-clean). Inert at 0 (geo
    // pristine; swap=0 fully pristine).
    //   geo_nonsep : 0 = off; 1 = on (non-separable rescore of the geo candidates).
    //   geo_ctarget: clearance cap in units of mean radius (0 => cap at 0 = contact; a
    //                positive value, e.g. 0.5, keeps g active on loose cages too).
    //   geo_nsk    : top-K multiplier over geo_B for the rescore set (default 4).
    int    geo_nonsep = 0;
    double geo_ctarget = 0.5;
    int    geo_nsk = 4;
    // [190] H26: step window for geo_nonsep==1 - the coupled rescore fires only for
    // nonsep_start <= current_step_ < nonsep_stop (stop < 0 = open). Default = always.
    long   geo_nonsep_start = 0;
    long   geo_nonsep_stop  = -1;
};

// Near-pair edge sink filled by get_forces_nd_3 (only on sweep steps, only when
// swap > 0): every pair with r < (Di+Dj)/2 + margin, per thread. The sweep needs
// exactly these pairs -- any pair farther than the largest admissible radius change
// contributes zero to every candidate dE.
struct SwapEdge { std::uint32_t i, j; double d2; };
struct SwapEdgeSink {
    // --- exact mode: near-pair edges ---
    bool want_edges = false;
    double margin = 0.0;
    // Per-thread buffers, preallocated to `capacity` edges so the force loop can
    // append branch-free (write, then advance by the predicate). count[t] = used.
    std::vector<std::vector<SwapEdge>> per_thread;
    std::vector<std::size_t> count;
    // --- proxy mode: per-particle sum of absolute overlaps (contact - separation) ---
    bool want_gsum = false;
    std::vector<double> gsum_thread;   // nthreads * gstride; only each thread's touched
                                       // range is dirty between calls (see g_lo/g_hi)
    std::vector<double> gsum;          // reduced, size N
    std::vector<std::size_t> g_lo, g_hi;   // per-thread touched particle range of the LAST g write
    std::size_t gstride = 0;
    void reset(int nthreads, std::size_t capacity, std::size_t N) {
        if (want_edges) {
            if (per_thread.size() != static_cast<std::size_t>(nthreads)) per_thread.resize(nthreads);
            count.assign(nthreads, 0);
            for (auto& v : per_thread) if (v.size() < capacity + 8) v.resize(capacity + 8);
        }
        if (want_gsum) {
            const std::size_t stride = ((N + 7) / 8) * 8;
            const std::size_t need = static_cast<std::size_t>(nthreads) * stride;
            if (gstride != stride || gsum_thread.size() != need ||
                g_lo.size() != static_cast<std::size_t>(nthreads)) {
                gstride = stride;
                gsum_thread.assign(need, 0.0);          // full zero on (re)allocation
                g_lo.assign(nthreads, 0); g_hi.assign(nthreads, 0);
            } else {
                for (int t = 0; t < nthreads; ++t) {    // zero only what the last call dirtied
                    double* g = gsum_thread.data() + static_cast<std::size_t>(t) * gstride;
                    for (std::size_t i = g_lo[t]; i < g_hi[t]; ++i) g[i] = 0.0;
                }
            }
            gsum.assign(N, 0.0);
        }
    }
};

struct SwapSweepStats {
    bool starved = false;       // ctrl: quota unmet with the window at its floor
    double w_used = 0.0;        // ctrl: window in effect during this sweep
    std::uint64_t candidates = 0;
    std::uint64_t evaluated = 0;   // candidates that passed the O(1) prefilter
    std::uint64_t accepted = 0;
    double dE_sum = 0.0;        // total objective change (negative = lowered)
    double e_scale = 0.0;       // mean per-contact energy before the sweep
    std::size_t contacts = 0;
    // rescue channel (this sweep): decision cross-table of proxy vs exact
    std::uint64_t rescue_candidates = 0;
    std::uint64_t rescue_accepted = 0;      // exact says yes (move applied)
    std::uint64_t rescue_false_rej = 0;     // proxy rejects, exact accepts
    std::uint64_t rescue_false_acc = 0;     // proxy accepts, exact rejects
    // 3-cycles (this pass)
    std::uint64_t cycle3_candidates = 0;
    std::uint64_t cycle3_accepted = 0;
    std::uint64_t cycle3_blocked = 0;       // accepted cycles with ALL pairwise swaps uphill
};

class SwapEngine {
public:
    SwapEngine(const SwapParams& p, std::uint64_t seed);

    // Call whenever the index<->particle mapping changed (spatial reorder) or once
    // before the first sweep. Rebuilds the diameter-rank tables from D0.
    void rebuild_ranks(const std::vector<double>& D0);

    // Margin the force loop must use when collecting near pairs for the next
    // sweep: cap * (largest radius change any admissible candidate can make).
    // Rebuilds the rank tables if they were invalidated.
    double edge_margin(const std::vector<double>& D, const std::vector<double>& D0);

    // One sweep at fixed positions over the near-pair edges collected by the force
    // loop for the current x/D. Mutates D and D0 for accepted swaps (both arrays
    // exchange entries i<->j), appends every touched particle index to `touched`,
    // and returns statistics.
    SwapSweepStats sweep(std::size_t N,
                         std::vector<double>& D,
                         std::vector<double>& D0,
                         const SwapEdgeSink& edges,
                         const std::vector<std::size_t>& z,   // per-particle contact count (force loop)
                         std::vector<std::uint32_t>& touched);

    // Checkpoint support [log 043]: dump/restore the engine's full mutable state
    // (RNG, adaptive window, rank + window tables AS BUILT, counters, totals) to
    // an open binary stream. The tables are serialized rather than rebuilt: the
    // window table lags w_cur_ by up to 2% by design, so a rebuild-on-load would
    // diverge from an unbroken run. Returns false on short write/read.
    bool save_state(std::FILE* fp) const;
    bool load_state(std::FILE* fp);

    // Surgical injection (ME0-ii [log 054]): apply the trade (i, j) with
    // bookkeeping identical to a sweep acceptance (ranks, counts, touched,
    // event log with relief = NaN). Returns false on invalid indices.
    std::size_t event_width() const { return p_.event_log_geom ? 12u : 8u; }
    std::uint64_t geo_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                           const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched);
    // H2 [107/117]: 3-cycle pass, run AFTER geo_pass when p_.geo_cycle != 0.
    std::uint64_t geo_cycle_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                 const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched);
    // H34 [194]: 4-cycle pass, run AFTER geo_cycle_pass when p_.geo_cycle4 != 0.
    std::uint64_t geo_cycle4_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                  const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched);
    // H7 [121/122]: cage-capacity pass, run AFTER geo_pass when p_.geo_capk != 0.
    // H8 [129/130/131]: overlap-energy pairwise pass, run AFTER geo_pass when p_.geo_cage != 0.
    // x/box/Ndim are used ONLY when geo_cage_move != 0 (the diagnostic position move);
    // pass x=nullptr for the pure-permutation (C5-valid) form.
    std::uint64_t geo_cage_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched,
                                double* x = nullptr, const double* box = nullptr, int Ndim = 3);
    // H12 [146/148]: bottleneck-relief pass (reopen-loop "final squeeze").
    std::uint64_t geo_bottleneck_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                      const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched);
    std::uint64_t geo_capk_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched);
    std::uint64_t geo_trades_ = 0, geo_passes_ = 0;
    double geo_T_cur_ = -1.0;       // [135] live annealed temperature (lazily init to geo_T)
    std::uint64_t geo_cycle_trades_ = 0;
    std::uint64_t geo_cycle4_trades_ = 0;   // [194] H34
    std::uint64_t geo_capk_trades_ = 0;
    std::uint64_t geo_cage_trades_ = 0;
    std::uint64_t geo_bottleneck_trades_ = 0;
    // [172] H19: set by the driver before a geo_pass; when geo_nonsep==2 the non-separable
    // rescore fires ONLY while this is true (i.e. only at the terminal reopen), so c_min
    // does the compression and the coupled search probes only the converged terminal.
    bool geo_reopen_active_ = false;
    bool inject(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                std::vector<std::uint32_t>& touched,
                std::uint32_t i, std::uint32_t j);

    const SwapParams& params() const { return p_; }
    // accepted swaps per particle index over the run (diagnostic; follows the index
    // through spatial reorders because the caller invalidates ranks -> we re-key by
    // permutation? no: counts are per current index; the caller permutes them)
    std::uint64_t swap_count(std::size_t k) const { return k < swap_count_.size() ? swap_count_[k] : 0; }
    std::vector<std::uint64_t>& swap_counts() { return swap_counts_ref(); }
    bool ranks_valid() const { return ranks_valid_; }
    void invalidate_ranks() { ranks_valid_ = false; }

private:
    SwapParams p_;
    std::uint64_t rng_;
    bool ranks_valid_ = false;
    std::vector<std::uint32_t> rank_of_;   // particle -> rank (ascending D0)
    std::vector<std::uint32_t> at_rank_;   // rank -> particle
    // mode 5: per-rank admissible partner-rank window for |D_j/D_i - 1| <= window_rel.
    // The diameter multiset is fixed, so this table depends only on the sorted D0.
    std::vector<std::uint32_t> win_lo_, win_hi_;   // inclusive
    // rank-space working arrays for the proxy sweep (streaming access): value at
    // rank r belongs to particle at_rank_[r]
    std::vector<double> rank_D0_;          // sorted D0 (fixed multiset; rebuilt with ranks)
    std::size_t wmax_cached_ = 0;          // max window span in ranks (0 = recompute)
    double w_cur_ = -1.0;                  // ctrl: adaptive window state (-1 = uninitialised)
    void rebuild_window_table(double w);   // two-pointer over rank_D0_
    void build_window_table_for(double w, std::vector<std::uint32_t>& lo,
                                std::vector<std::uint32_t>& hi) const;
    // rescue channel: its own (wider) window table over the same sorted D0; the
    // table depends only on the sorted multiset + rescue_w, so it is built lazily
    // once per rank rebuild
    std::vector<std::uint32_t> rwin_lo_, rwin_hi_;
    bool rwin_built_ = false;
    void rescue_pass(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                     const SwapEdgeSink& edges, std::vector<std::uint32_t>& touched,
                     SwapSweepStats& st);
    std::vector<double> rank_g_;
    std::vector<std::uint32_t> rank_z_, rank_ppc_;
    // scratch (persist across sweeps to avoid reallocation)
    std::vector<std::uint32_t> deg_, adj_off_, adj_idx_, fill_;
    std::vector<double> adj_r_;
    std::vector<double> rad_;
    std::vector<double> cage_, fact_;      // ME0-iv per-particle geometry (event log only)
    void geom_from_edges(std::size_t N, const std::vector<double>& D,
                         const SwapEdgeSink& edges, const std::vector<std::size_t>& z);
    std::vector<double> gsum_;             // per-particle sum of current overlaps (prefilter)
    std::vector<std::uint32_t> zc_;        // per-particle current contact count (prefilter)
    std::vector<std::uint32_t> per_particle_count_;
    std::vector<std::uint64_t> swap_count_;   // per current index
    std::vector<std::uint64_t>& swap_counts_ref() { return swap_count_; }
public:
    // profiling (RCP_PROFILE): accumulated seconds per phase and counters
    double t_build_ = 0.0, t_cand_ = 0.0;
    std::uint64_t n_edges_total_ = 0, n_eval_total_ = 0, n_cand_total_ = 0, n_sweeps_ = 0;
    // per-trade event log (ME0-i): caller sets current_step_ before each sweep
    std::uint64_t current_step_ = 0;
    std::vector<double> event_buf_;
    bool event_log_full_ = false;
    // rescue channel run totals (mirrored into PackingResult)
    std::uint64_t rescue_cand_total_ = 0, rescue_acc_total_ = 0;
    std::uint64_t rescue_frej_total_ = 0, rescue_facc_total_ = 0;
    std::uint64_t cycle3_cand_total_ = 0, cycle3_acc_total_ = 0, cycle3_blocked_total_ = 0;
private:

    std::uint64_t next_u64();
    double next_unit();
    double delta_E(std::uint32_t i, std::uint32_t j) const;
    // exact dE of the rotation i<-R_j, j<-R_k, k<-R_i (intra-trio pairs included)
    double delta_E3(std::uint32_t i, std::uint32_t j, std::uint32_t k) const;
    SwapSweepStats sweep_proxy(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                               const double* gsum_in, const std::vector<std::size_t>& z,
                               std::vector<std::uint32_t>& touched);
    SwapSweepStats sweep_proxy_rank(std::size_t N, std::vector<double>& D, std::vector<double>& D0,
                                    const double* gsum_in, const std::vector<std::size_t>& z,
                                    std::vector<std::uint32_t>& touched);
};

}  // namespace rcpgenerator
