#include "rcpgenerator/rcp_generator.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <omp.h>

// libstdc++ parallel-mode <parallel/algorithm> provides __gnu_parallel::sort,
// an OpenMP-backed quicksort/mergesort. Used to parallelize the index sort by
// x-coordinate inside get_pairs_nd_3 (previously serial and a measured ~15-40%
// of get_pairs cost at N=50K-100K). Falls back to std::sort if unavailable.
#if defined(__GNUC__) && __has_include(<parallel/algorithm>)
#  include <parallel/algorithm>
#  define RCP_HAVE_GNU_PARALLEL_SORT 1
#else
#  define RCP_HAVE_GNU_PARALLEL_SORT 0
#endif

// Cycle 15a: cycle-counter instrumentation (rdtsc). Lets us measure
// cycles-per-candidate inside the inner force loop and get_pairs walk to
// quantify L1/L2/L3 contention. Toggled via env var RCP_CYCLE_PROBE=1.
#include <x86intrin.h>
#include <cstdlib>

static std::atomic<std::uint64_t> g_floop_cycles{0};
static std::atomic<std::uint64_t> g_floop_candidates{0};
static std::atomic<std::uint64_t> g_floop_overlaps{0};
// Cycle 15c: split forces.loop cycles into test-path (always-executed)
// and overlap-path (sqrt + buffer writes, ~3.7% rate).
static std::atomic<std::uint64_t> g_floop_test_cycles{0};
static std::atomic<std::uint64_t> g_floop_overlap_cycles{0};
static std::atomic<std::uint64_t> g_gpairs_cycles{0};
static std::atomic<std::uint64_t> g_gpairs_candidates{0};
static std::atomic<std::uint64_t> g_gpairs_accepted{0};
// Cycle 15c: split get_pairs cycles into cell-traversal, distance-test,
// and accepted-write phases.
static std::atomic<std::uint64_t> g_gpairs_cell_cycles{0};
static std::atomic<std::uint64_t> g_gpairs_dist_cycles{0};
static std::atomic<std::uint64_t> g_gpairs_write_cycles{0};
// Cycle 15i: measure the per-call sort cost in get_pairs.
static std::atomic<std::uint64_t> g_gpairs_sort_cycles{0};
static std::atomic<std::uint64_t> g_gpairs_sort_calls{0};
// Cycle 15c: per-thread wall-time for get_pairs, to expose thread imbalance.
static constexpr int RCP_PROBE_MAX_THREADS = 16;
static std::atomic<std::uint64_t> g_gpairs_thread_cycles[RCP_PROBE_MAX_THREADS] = {};

// kdtree Cycle 4: per-phase RDTSC probes for get_pairs_kdtree, mirroring the
// g_gpairs_* breakdown used to drive sort+walk optimization. Each probe is
// accumulated thread-locally inside the parallel region and folded into the
// atomic at function exit.
static std::atomic<std::uint64_t> g_kd_total_cycles{0};   // wall covering get_pairs_kdtree (master thread)
static std::atomic<std::uint64_t> g_kd_setup_cycles{0};   // Dmax/median/shell + inv_box + j_locks init
static std::atomic<std::uint64_t> g_kd_build_cycles{0};   // tree.build
static std::atomic<std::uint64_t> g_kd_image_cycles{0};   // compute_image_set + xq construction per image
static std::atomic<std::uint64_t> g_kd_query_cycles{0};   // tree.range_query
static std::atomic<std::uint64_t> g_kd_mic_cycles{0};     // per-candidate MIC d^2 test + r_c^2 compare
static std::atomic<std::uint64_t> g_kd_write_cycles{0};   // accept-write + dup check + locks
static std::atomic<std::uint64_t> g_kd_sort_cycles{0};    // per-particle final sort
static std::atomic<std::uint64_t> g_kd_calls{0};
static std::atomic<std::uint64_t> g_kd_particles_queried{0};
static std::atomic<std::uint64_t> g_kd_candidates{0};
static std::atomic<std::uint64_t> g_kd_accepted{0};
static std::atomic<std::uint64_t> g_kd_thread_cycles[RCP_PROBE_MAX_THREADS] = {};

// Cycle 7 measurement: histograms of refresh count + leaves touched per call.
// Buckets are log-spaced powers of 2 with overflow into last bucket.
//   bin 0: refresh=0
//   bin k (k>=1): 2^(k-1) <= refresh < 2^k, up to bin 18 (262144+)
static constexpr int RCP_KD_HIST_BINS = 19;
static std::atomic<std::uint64_t> g_kd_refresh_hist[RCP_KD_HIST_BINS] = {};
static std::atomic<std::uint64_t> g_kd_leaves_touched_hist[RCP_KD_HIST_BINS] = {};
// Per-regime breakdown: 0=early(steps 0-2499), 1=mid(2500-9999), 2=late(10000+)
static constexpr int RCP_KD_REGIME_COUNT = 3;
static std::atomic<std::uint64_t> g_kd_regime_calls[RCP_KD_REGIME_COUNT] = {};
static std::atomic<std::uint64_t> g_kd_regime_refresh_sum[RCP_KD_REGIME_COUNT] = {};
static std::atomic<std::uint64_t> g_kd_regime_leaves_touched_sum[RCP_KD_REGIME_COUNT] = {};
static std::atomic<std::uint64_t> g_kd_regime_build_cyc[RCP_KD_REGIME_COUNT] = {};
static std::atomic<std::uint64_t> g_kd_regime_query_cyc[RCP_KD_REGIME_COUNT] = {};
static int rcp_kd_regime_for_step(std::uint64_t step) {
    if (step < 2500)  return 0;
    if (step < 10000) return 1;
    return 2;
}
static int rcp_kd_log2_bin(std::uint64_t v) {
    if (v == 0) return 0;
    int b = 1;
    while (v > 1 && b < RCP_KD_HIST_BINS - 1) { v >>= 1; ++b; }
    return b;
}

static bool rcp_cycle_probe_enabled() {
    static int v = -1;
    if (v == -1) {
        const char* s = std::getenv("RCP_CYCLE_PROBE");
        v = (s && s[0] && s[0] != '0') ? 1 : 0;
    }
    return v != 0;
}
// Cycle 15c MIC-cost probe: number of additional MIC ops per dim to inject
// into the inner force loop. Their results sink into an atomic global so the
// compiler can't elide them. Slope of cyc/cand vs this count = per-MIC cost.
static int rcp_extra_mic_count() {
    static int v = -1;
    if (v == -1) {
        const char* s = std::getenv("RCP_EXTRA_MIC");
        v = s ? std::atoi(s) : 0;
        if (v < 0) v = 0;
        if (v > 16) v = 16;
    }
    return v;
}
// Cycle 25/26: 16-byte packed entry for get_pairs cell-walk. Holds j,
// the float-precision D[j]*0.5 (used in r_c computation), and x[j*Ndim+0].
// Reading packed_sort[pos] gets all three in one cache line — replaces
// separate sort_idx[pos] + D[j] + x[j*Ndim] loads that dominated
// cell-walk per-cand cost.
//
// D_half_f is float (~7 decimal digits) — losing ~1 ulp of precision on
// r_c. Borderline shell-edge pairs may classify differently within ~1e-7;
// physically negligible for size distributions where D ~ O(1).
struct PackedSortEntry {
    std::uint32_t j;
    float         D_half_f;   // D[j] * 0.5 as float
    double        x0;
};
static_assert(sizeof(PackedSortEntry) == 16, "PackedSortEntry must be 16 bytes");
// ============================================================================
// KDTREE PROJECT — Cycle 1 scaffolding: file-scope infrastructure
// ----------------------------------------------------------------------------
// Backend selector (via RCP_PAIRS_BACKEND env var), shadow-mode counters,
// and reset/dump helpers. These live at file scope so they can be called from
// rcp_cycle_probe_reset() (also file-scope). The actual get_pairs_kdtree and
// get_pairs_dispatch live inside namespace rcpgenerator alongside
// get_pairs_nd_3 (see further down).
// ----------------------------------------------------------------------------

enum class PairsBackend { SortWalk, KdTree, Shadow };

// Cycle 10 H43: select binary KdTree (default) vs OctTree for the
// neighbor-finding path under PairsBackend::KdTree. Read once.
static bool rcp_use_octree() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_KDTREE_VARIANT");
        return s && std::strcmp(s, "octree") == 0;
    }();
    return v;
}
// Cycle 13 b6.A Phase 2: read K_BATCH so the kdtree patch-path inflation
// can scale with cumulative-drift bound (K × Dmin/4 instead of fixed Dmin/2).
static int rcp_k_batch() {
    static const int v = []() {
        const char* s = std::getenv("RCP_K_BATCH");
        if (!s || !s[0]) return 1;
        int x = std::atoi(s);
        return (x >= 1 && x <= 16) ? x : 1;
    }();
    return v;
}
static PairsBackend rcp_pairs_backend() {
    static const PairsBackend v = []() {
        const char* s = std::getenv("RCP_PAIRS_BACKEND");
        if (!s || !*s) return PairsBackend::SortWalk;
        if (std::strcmp(s, "kdtree") == 0) return PairsBackend::KdTree;
        if (std::strcmp(s, "shadow") == 0) return PairsBackend::Shadow;
        return PairsBackend::SortWalk;
    }();
    return v;
}

static std::atomic<std::uint64_t> g_shadow_calls{0};
static std::atomic<std::uint64_t> g_shadow_particles_diff{0};
static std::atomic<std::uint64_t> g_shadow_pairs_missing_in_kdtree{0};
static std::atomic<std::uint64_t> g_shadow_pairs_extra_in_kdtree{0};
static std::atomic<std::uint64_t> g_shadow_first_diff_step{UINT64_MAX};
static std::atomic<std::uint64_t> g_shadow_first_diff_particle{UINT64_MAX};
static thread_local std::uint64_t g_shadow_current_step = 0;

static void rcp_shadow_set_step(std::uint64_t step) {
    g_shadow_current_step = step;
}

static void rcp_shadow_reset() {
    g_shadow_calls = 0;
    g_shadow_particles_diff = 0;
    g_shadow_pairs_missing_in_kdtree = 0;
    g_shadow_pairs_extra_in_kdtree = 0;
    g_shadow_first_diff_step = UINT64_MAX;
    g_shadow_first_diff_particle = UINT64_MAX;
}

static void rcp_shadow_dump() {
    if (rcp_pairs_backend() != PairsBackend::Shadow) return;
    const std::uint64_t calls = g_shadow_calls.load();
    if (calls == 0) return;
    std::cerr << "\n=== Shadow-mode neighbor-list comparison ===\n";
    std::cerr << "  total get_pairs calls compared: " << calls << "\n";
    std::cerr << "  total particles with diff:      " << g_shadow_particles_diff.load() << "\n";
    std::cerr << "  pairs missing in kdtree:        " << g_shadow_pairs_missing_in_kdtree.load() << "\n";
    std::cerr << "  pairs extra in kdtree:          " << g_shadow_pairs_extra_in_kdtree.load() << "\n";
    const std::uint64_t fds = g_shadow_first_diff_step.load();
    if (fds != UINT64_MAX) {
        std::cerr << "  first diff step:                " << fds << "\n";
        std::cerr << "  first diff particle:            "
                  << g_shadow_first_diff_particle.load() << "\n";
    } else {
        std::cerr << "  STATUS: no divergence detected (set-equal across all calls)\n";
    }
}
// ============================================================================
// end kdtree file-scope scaffolding
// ============================================================================

static void rcp_cycle_probe_reset() {
    g_floop_cycles = 0;
    g_floop_candidates = 0;
    g_floop_overlaps = 0;
    g_floop_test_cycles = 0;
    g_floop_overlap_cycles = 0;
    g_gpairs_cycles = 0;
    g_gpairs_candidates = 0;
    g_gpairs_accepted = 0;
    g_gpairs_cell_cycles = 0;
    g_gpairs_dist_cycles = 0;
    g_gpairs_write_cycles = 0;
    g_gpairs_sort_cycles = 0;
    g_gpairs_sort_calls = 0;
    for (int t = 0; t < RCP_PROBE_MAX_THREADS; ++t) {
        g_gpairs_thread_cycles[t] = 0;
    }
    g_kd_total_cycles = 0;
    g_kd_setup_cycles = 0;
    g_kd_build_cycles = 0;
    g_kd_image_cycles = 0;
    g_kd_query_cycles = 0;
    g_kd_mic_cycles = 0;
    g_kd_write_cycles = 0;
    g_kd_sort_cycles = 0;
    g_kd_calls = 0;
    g_kd_particles_queried = 0;
    g_kd_candidates = 0;
    g_kd_accepted = 0;
    for (int t = 0; t < RCP_PROBE_MAX_THREADS; ++t) {
        g_kd_thread_cycles[t] = 0;
    }
    for (int b = 0; b < RCP_KD_HIST_BINS; ++b) {
        g_kd_refresh_hist[b] = 0;
        g_kd_leaves_touched_hist[b] = 0;
    }
    for (int r = 0; r < RCP_KD_REGIME_COUNT; ++r) {
        g_kd_regime_calls[r] = 0;
        g_kd_regime_refresh_sum[r] = 0;
        g_kd_regime_leaves_touched_sum[r] = 0;
        g_kd_regime_build_cyc[r] = 0;
        g_kd_regime_query_cyc[r] = 0;
    }
    rcp_shadow_reset();
}
static void rcp_cycle_probe_dump() {
    if (!rcp_cycle_probe_enabled()) return;
    std::uint64_t fc = g_floop_cycles.load();
    std::uint64_t fn = g_floop_candidates.load();
    std::uint64_t fo = g_floop_overlaps.load();
    std::uint64_t ft = g_floop_test_cycles.load();
    std::uint64_t fv = g_floop_overlap_cycles.load();
    std::uint64_t gc = g_gpairs_cycles.load();
    std::uint64_t gn = g_gpairs_candidates.load();
    std::uint64_t ga = g_gpairs_accepted.load();
    std::uint64_t gcell = g_gpairs_cell_cycles.load();
    std::uint64_t gdist = g_gpairs_dist_cycles.load();
    std::uint64_t gwrite = g_gpairs_write_cycles.load();
    std::cerr << "\n=== RDTSC cycle-probe (Cycle 15c) ===\n";
    if (fn > 0) {
        std::cerr << "  forces.loop:  " << fc << " cyc  "
                  << fn << " cand  "
                  << fo << " ovlp  ="
                  << "  " << (double)fc / fn << " cyc/cand"
                  << "  " << (double)fc / std::max<std::uint64_t>(fo,1) << " cyc/ovlp\n";
        std::cerr << "    test-path:    " << ft << " cyc  "
                  << (double)ft / fn << " cyc/cand"
                  << "  (" << 100.0 * ft / std::max<std::uint64_t>(fc,1) << "%)\n";
        std::cerr << "    overlap-path: " << fv << " cyc  "
                  << (double)fv / std::max<std::uint64_t>(fo,1) << " cyc/ovlp"
                  << "  amortized " << (double)fv / fn << " cyc/cand"
                  << "  (" << 100.0 * fv / std::max<std::uint64_t>(fc,1) << "%)\n";
    }
    if (gn > 0) {
        std::cerr << "  get_pairs:    " << gc << " cyc  "
                  << gn << " cand  "
                  << ga << " acc  ="
                  << "  " << (double)gc / gn << " cyc/cand\n";
        if (gcell + gdist + gwrite > 0) {
            std::cerr << "    cell-walk:    " << gcell << " cyc  ("
                      << 100.0 * gcell / std::max<std::uint64_t>(gc,1) << "%)\n";
            std::cerr << "    dist-test:    " << gdist << " cyc  ("
                      << 100.0 * gdist / std::max<std::uint64_t>(gc,1) << "%)\n";
            std::cerr << "    accept-write: " << gwrite << " cyc  ("
                      << 100.0 * gwrite / std::max<std::uint64_t>(gc,1) << "%)\n";
        }
        {
            std::uint64_t gsc = g_gpairs_sort_cycles.load();
            std::uint64_t gscn = g_gpairs_sort_calls.load();
            if (gscn > 0) {
                std::cerr << "    sort:         " << gsc << " cyc  "
                          << gscn << " calls  avg "
                          << (double)gsc / gscn << " cyc/call  ("
                          << 100.0 * gsc / std::max<std::uint64_t>(gc,1) << "% of get_pairs)\n";
            }
        }
        std::cerr << "    per-thread cyc:";
        for (int t = 0; t < RCP_PROBE_MAX_THREADS; ++t) {
            std::uint64_t tc = g_gpairs_thread_cycles[t].load();
            if (tc > 0) std::cerr << " t" << t << "=" << tc;
        }
        std::cerr << "\n";
    }
    // kdtree Cycle 4: per-phase breakdown of get_pairs_kdtree.
    std::uint64_t kdcalls = g_kd_calls.load();
    if (kdcalls > 0) {
        std::uint64_t ktotal = g_kd_total_cycles.load();
        std::uint64_t ksetup = g_kd_setup_cycles.load();
        std::uint64_t kbuild = g_kd_build_cycles.load();
        std::uint64_t kimage = g_kd_image_cycles.load();
        std::uint64_t kquery = g_kd_query_cycles.load();
        std::uint64_t kmic   = g_kd_mic_cycles.load();
        std::uint64_t kwrite = g_kd_write_cycles.load();
        std::uint64_t ksort  = g_kd_sort_cycles.load();
        std::uint64_t kparts = g_kd_particles_queried.load();
        std::uint64_t kcand  = g_kd_candidates.load();
        std::uint64_t kacc   = g_kd_accepted.load();
        std::cerr << "  kdtree:       " << ktotal << " cyc (master-wall sum, "
                  << kdcalls << " calls)\n";
        std::cerr << "    setup:        " << ksetup << " cyc  ("
                  << 100.0 * ksetup / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    build:        " << kbuild << " cyc  ("
                  << 100.0 * kbuild / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    image-enum:   " << kimage << " cyc  ("
                  << 100.0 * kimage / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    range-query:  " << kquery << " cyc  ("
                  << 100.0 * kquery / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    mic-test:     " << kmic   << " cyc  ("
                  << 100.0 * kmic   / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    accept-write: " << kwrite << " cyc  ("
                  << 100.0 * kwrite / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        std::cerr << "    sort:         " << ksort  << " cyc  ("
                  << 100.0 * ksort  / std::max<std::uint64_t>(ktotal,1) << "%)\n";
        if (kparts > 0) {
            std::cerr << "    queried particles=" << kparts
                      << "  candidates=" << kcand
                      << "  accepted="    << kacc
                      << "  cand/part="   << (double)kcand / kparts
                      << "  acc/cand="    << (double)kacc  / std::max<std::uint64_t>(kcand,1) << "\n";
        }
        std::cerr << "    per-thread cyc:";
        for (int t = 0; t < RCP_PROBE_MAX_THREADS; ++t) {
            std::uint64_t tc = g_kd_thread_cycles[t].load();
            if (tc > 0) std::cerr << " t" << t << "=" << tc;
        }
        std::cerr << "\n";
        // Cycle 7 measurement dump: refresh / leaves-touched histogram +
        // per-regime breakdown to validate the model assumptions.
        std::cerr << "    refresh-count histogram (calls by refresh count, log2 bins):\n";
        std::cerr << "      ";
        for (int b = 0; b < RCP_KD_HIST_BINS; ++b) {
            std::uint64_t cnt = g_kd_refresh_hist[b].load();
            if (cnt == 0) continue;
            std::uint64_t lo = (b == 0) ? 0 : (1ull << (b-1));
            std::uint64_t hi = (b == 0) ? 0 : ((1ull << b) - 1);
            std::cerr << "[" << lo << "-" << hi << "]=" << cnt << " ";
        }
        std::cerr << "\n";
        std::cerr << "    leaves-touched histogram (calls by # distinct leaves with refresh):\n";
        std::cerr << "      ";
        for (int b = 0; b < RCP_KD_HIST_BINS; ++b) {
            std::uint64_t cnt = g_kd_leaves_touched_hist[b].load();
            if (cnt == 0) continue;
            std::uint64_t lo = (b == 0) ? 0 : (1ull << (b-1));
            std::uint64_t hi = (b == 0) ? 0 : ((1ull << b) - 1);
            std::cerr << "[" << lo << "-" << hi << "]=" << cnt << " ";
        }
        std::cerr << "\n";
        std::cerr << "    per-regime (early=0-2499, mid=2500-9999, late=10000+):\n";
        const char* regime_name[] = {"early", "mid  ", "late "};
        for (int r = 0; r < RCP_KD_REGIME_COUNT; ++r) {
            std::uint64_t rc = g_kd_regime_calls[r].load();
            if (rc == 0) continue;
            std::uint64_t rs = g_kd_regime_refresh_sum[r].load();
            std::uint64_t lt = g_kd_regime_leaves_touched_sum[r].load();
            std::uint64_t bc = g_kd_regime_build_cyc[r].load();
            std::uint64_t qc = g_kd_regime_query_cyc[r].load();
            std::cerr << "      " << regime_name[r]
                      << " calls=" << rc
                      << " avg_refresh=" << (double)rs/rc
                      << " avg_leaves_touched=" << (double)lt/rc
                      << " build_cyc/call=" << bc/rc
                      << " query_cyc/call=" << qc/rc << "\n";
        }
    }
}

namespace rcpgenerator {

namespace {

static constexpr double ALPHA_MAX = 0.0025;
static constexpr double BETA1 = 0.9;
static constexpr double BETA2 = 0.999;
static constexpr double EPSILON = 1e-8;
static constexpr std::size_t N_STEPS = 60000;
static constexpr double DT = 0.1;
static const std::string METHOD = "ADAM";
static const std::vector<std::uint32_t> MAX_NEIGHBORS = {300, 750, 5500, 5500};
static constexpr double DELTA_PHI0 = 1.5 * 1.5e-3;

// ===== profile instrumentation (added for kernel-timing audit) =====
struct KernelTimer {
    std::chrono::steady_clock::time_point start;
    double total_seconds = 0.0;
    std::size_t call_count = 0;
    void begin() { start = std::chrono::steady_clock::now(); }
    void end() {
        auto dt = std::chrono::steady_clock::now() - start;
        total_seconds += std::chrono::duration<double>(dt).count();
        ++call_count;
    }
    void reset() { total_seconds = 0.0; call_count = 0; }
};

static KernelTimer g_t_pairs;
static KernelTimer g_t_forces;
static KernelTimer g_t_adam;
static KernelTimer g_t_position;
static KernelTimer g_t_refresh;
static KernelTimer g_t_bookkeep;  // everything else inside the main loop

// Sub-timers inside get_forces_nd_3 — for Cycle 3 scaling diagnosis.
// Total of these < g_t_forces (begin/end overhead).
static KernelTimer g_t_forces_alloc;     // F_local_flat / z_local_flat assign
static KernelTimer g_t_forces_xflat;     // shadow refresh of x_flat
static KernelTimer g_t_forces_pflat;     // shadow refresh of pairs_flat
static KernelTimer g_t_forces_loop;      // main parallel force loop + reductions
static KernelTimer g_t_forces_reduce;    // per-thread F/z merge into shared

// Cycle 15n: sub-timers inside bookkeep to find which sub-region is heavy.
// These overlap with g_t_bookkeep (count subsets of bookkeep time).
static KernelTimer g_t_bk_dupdate;      // D[i] = D0[i] * kappa loop
static KernelTimer g_t_bk_meanforce;    // compute_mean_force + F_history
static KernelTimer g_t_bk_alpha;        // alpha trend logic + history scans
static KernelTimer g_t_bk_phistory;     // mu_flag phi_history scans
static KernelTimer g_t_bk_misc;         // catch-all per-step scalar work

// Version counters for pflat shadow refresh skip (Cycle 3b). get_pairs_nd_3
// increments g_pairs_version after every pair-list update. get_forces_nd_3
// compares this against g_pflat_version and refreshes pairs_flat only when
// they differ. Reset at the top of run_packing_observed so versions
// don't leak across runs.
static std::uint64_t g_pairs_version = 0;
static std::uint64_t g_pflat_version = 0;

static void print_timing_summary(double total_wall, std::size_t steps) {
    std::cerr << "\n=== Per-kernel timing breakdown ===\n";
    auto show = [&](const char* name, const KernelTimer& t) {
        double pct = total_wall > 0 ? 100.0 * t.total_seconds / total_wall : 0.0;
        double avg_ms = t.call_count ? (t.total_seconds / t.call_count) * 1000.0 : 0.0;
        std::cerr << "  " << name
                  << "  total=" << t.total_seconds << "s"
                  << "  (" << pct << "%)"
                  << "  calls=" << t.call_count
                  << "  avg=" << avg_ms << "ms\n";
    };
    show("get_pairs_nd_3   ", g_t_pairs);
    show("get_forces_nd_3  ", g_t_forces);
    show("  forces.alloc   ", g_t_forces_alloc);
    show("  forces.xflat   ", g_t_forces_xflat);
    show("  forces.pflat   ", g_t_forces_pflat);
    show("  forces.loop    ", g_t_forces_loop);
    show("  forces.reduce  ", g_t_forces_reduce);
    show("adam_update      ", g_t_adam);
    show("position+pbc     ", g_t_position);
    show("refresh_check    ", g_t_refresh);
    show("bookkeep         ", g_t_bookkeep);
    show("  bk.D-update    ", g_t_bk_dupdate);
    show("  bk.meanforce   ", g_t_bk_meanforce);
    show("  bk.alpha       ", g_t_bk_alpha);
    show("  bk.phistory    ", g_t_bk_phistory);
    show("  bk.misc        ", g_t_bk_misc);
    double accounted = g_t_pairs.total_seconds + g_t_forces.total_seconds
                     + g_t_adam.total_seconds + g_t_position.total_seconds
                     + g_t_refresh.total_seconds + g_t_bookkeep.total_seconds;
    std::cerr << "  -- accounted: " << accounted << "s of " << total_wall << "s ("
              << (total_wall > 0 ? 100.0 * accounted / total_wall : 0.0) << "%)\n";
    std::cerr << "  -- steps: " << steps << "\n";
    // Cycle 14 H53: always-on force.loop rejection rate.
    {
        const std::uint64_t cand = g_floop_candidates.load();
        const std::uint64_t ovlp = g_floop_overlaps.load();
        if (cand > 0) {
            const double accept_rate = static_cast<double>(ovlp) / cand;
            std::cerr << "  [shell] floop_candidates=" << cand
                      << " overlaps=" << ovlp
                      << " accept_rate=" << accept_rate
                      << " reject_rate=" << (1.0 - accept_rate) << "\n";
        }
    }
}

}  // namespace

double delta_x(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& x_old,
    const std::vector<double>& D)
{
    std::size_t N = x.size();
    if (x_old.size() != N || D.size() != N) {
        throw std::invalid_argument("delta_x: mismatched sizes");
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double d2 = 0.0;
        for (std::size_t d = 0; d < x[i].size(); ++d) {
            double diff = x_old[i][d] - x[i][d];
            d2 += diff * diff;
        }
        double dist = std::sqrt(d2);
        sum += dist / D[i];
    }

    return sum / double(N);
}

double norm(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("norm: vector sizes differ");
    }
    double sum2 = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        sum2 += d * d;
    }
    return std::sqrt(sum2);
}

double mean(
    const std::vector<double>& v,
    std::size_t start,
    std::size_t end)
{
    double sum = 0.0;
    for (std::size_t k = start; k <= end; ++k) {
        sum += v[k];
    }
    return sum / double(end - start + 1);
}

double compute_mean_force(
    const double* F,
    std::size_t N,
    std::size_t Ndim,
    double Lc,
    const std::vector<std::size_t>& z)
{
    double sum_mag = 0.0;
    double sum_z = 0.0;

    std::uint32_t count = 0;

    // Cycle 15h: this used to be serial, costing ~2s of bookkeep wall at
    // T=4 N=100K because it runs every step over all N particles. Reduction
    // is associative-up-to-FP-noise, same as forces.reduce.
    #pragma omp parallel for reduction(+:sum_mag,count) schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t base = i * Ndim;
        double d2 = 0.0;
        for (std::size_t d = 0; d < Ndim; ++d) {
            double Fid = F[base + d];
            d2 += Fid * Fid;
        }
        if (d2 > 1.0E-16) {
            sum_mag += std::sqrt(d2);
            ++count;
        }
    }

    double mean_mag = sum_mag / double(count);

    (void)Lc;
    (void)z;
    return mean_mag / std::sqrt(double(Ndim));
}

double max_elem(const std::vector<std::vector<double>>& M) {
    if (M.empty()) {
        throw std::invalid_argument("maxElem: empty matrix");
    }
    double m = M[0][0];
    for (const auto& row : M) {
        if (row.empty()) continue;
        double row_max = *std::max_element(row.begin(), row.end());
        if (row_max > m) m = row_max;
    }
    return m;
}

double compute_max_min_dist(const std::vector<std::vector<double>>& min_dist) {
    double m = 0.0;
    for (auto& row : min_dist)
        for (double v : row)
            if (v > m) m = v;
    return m;
}

std::size_t compute_max_neighbors(const std::vector<std::vector<std::uint32_t>>& pairs) {
    std::size_t m = 0;
    for (auto& row : pairs)
        if (row[0] > m) m = row[0];
    return m;
}

// =====================================================================
// Cycle 16: run-state save/load (rung-0 endpoint snapshot for warm
// restart). Format is raw IEEE 754 binary — NO text encoding anywhere —
// so resume is bit-exact across save/load.
// File layout (all little-endian, native packing):
//   [4] magic "RCPS"
//   [4] version (uint32) = 1
//   [8] N (uint64)
//   [8] Ndim (uint64)
//   [8] step (uint64) — value of `step` at save time
//   [8] hist_used (uint64) — number of history entries actually filled
//   Scalar block (raw struct dump, packed):
//     double alpha, mu, kappa
//     double m_kappa, v_kappa, v_update_kappa, m_hat_kappa, v_hat_kappa
//     double t_adam            // ADAM bias-correction counter
//     double phi, phi_modifier
//     double phi_max_value
//     uint64_t phi_max_step
//     uint64_t mu_change, LastPhiUpdate, LastAlphaUpdate
//     int32_t  mu_flag
//     int32_t  _pad
//     double box[8]            // KD_MAX_NDIM upper bound
//   Arrays (raw doubles, in this exact order):
//     positions     : N*Ndim
//     diameters     : N    (DERIVED state; D = D0 * kappa each step)
//     D0            : N    (TRUE state; the reference array)
//     m_flat        : N*Ndim
//     v_flat        : N*Ndim
//     phi_history   : hist_used
//     F_history     : hist_used
//     U_history     : hist_used
// =====================================================================
struct RcpSavedState {
    std::uint64_t step = 0;
    std::uint64_t hist_used = 0;
    double alpha = 0, mu = 0, kappa = 0;
    double m_kappa = 0, v_kappa = 0, v_update_kappa = 0;
    double m_hat_kappa = 0, v_hat_kappa = 0;
    double t_adam = 0;
    double phi = 0, phi_modifier = 0;
    double phi_max_value = 0;
    std::uint64_t phi_max_step = 0;
    std::uint64_t mu_change = 0, LastPhiUpdate = 0, LastAlphaUpdate = 0;
    std::int32_t mu_flag = 1;
    double box[8] = {0};
    std::vector<double> positions, diameters, D0, m_flat, v_flat;
    std::vector<double> phi_history, F_history, U_history;
};

static bool rcp_save_state(const std::string& path,
                           std::size_t N, std::size_t Ndim,
                           const RcpSavedState& s) {
    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::cerr << "[rcp_save_state] fopen failed: " << path << "\n";
        return false;
    }
    auto W = [&](const void* p, std::size_t n) {
        return std::fwrite(p, 1, n, fp) == n;
    };
    bool ok = true;
    const char magic[5] = "RCPS";
    const std::uint32_t version = 2;
    const std::uint64_t N64 = N, Ndim64 = Ndim;
    ok &= W(magic, 4);
    ok &= W(&version, 4);
    ok &= W(&N64, 8);
    ok &= W(&Ndim64, 8);
    ok &= W(&s.step, 8);
    ok &= W(&s.hist_used, 8);
    // Scalar block
    ok &= W(&s.alpha, 8);
    ok &= W(&s.mu, 8);
    ok &= W(&s.kappa, 8);
    ok &= W(&s.m_kappa, 8);
    ok &= W(&s.v_kappa, 8);
    ok &= W(&s.v_update_kappa, 8);
    ok &= W(&s.m_hat_kappa, 8);
    ok &= W(&s.v_hat_kappa, 8);
    ok &= W(&s.t_adam, 8);
    ok &= W(&s.phi, 8);
    ok &= W(&s.phi_modifier, 8);
    ok &= W(&s.phi_max_value, 8);
    ok &= W(&s.phi_max_step, 8);
    ok &= W(&s.mu_change, 8);
    ok &= W(&s.LastPhiUpdate, 8);
    ok &= W(&s.LastAlphaUpdate, 8);
    std::int32_t mu_flag32 = s.mu_flag, pad32 = 0;
    ok &= W(&mu_flag32, 4);
    ok &= W(&pad32, 4);
    ok &= W(s.box, 8 * 8);
    // Arrays
    auto Wvec = [&](const std::vector<double>& v) {
        return v.empty() ? true : W(v.data(), v.size() * sizeof(double));
    };
    if (s.positions.size() != N * Ndim) ok = false;
    if (s.diameters.size() != N) ok = false;
    if (s.D0.size() != N) ok = false;
    if (s.m_flat.size() != N * Ndim) ok = false;
    if (s.v_flat.size() != N * Ndim) ok = false;
    if (s.phi_history.size() != s.hist_used) ok = false;
    if (s.F_history.size() != s.hist_used) ok = false;
    if (s.U_history.size() != s.hist_used) ok = false;
    ok &= Wvec(s.positions);
    ok &= Wvec(s.diameters);
    ok &= Wvec(s.D0);
    ok &= Wvec(s.m_flat);
    ok &= Wvec(s.v_flat);
    ok &= Wvec(s.phi_history);
    ok &= Wvec(s.F_history);
    ok &= Wvec(s.U_history);
    std::fclose(fp);
    if (!ok) {
        std::cerr << "[rcp_save_state] write failed: " << path << "\n";
    } else {
        std::cerr << "[rcp_save_state] saved step=" << s.step
                  << " hist_used=" << s.hist_used
                  << " path=" << path << "\n";
    }
    return ok;
}

static bool rcp_load_state(const std::string& path,
                           std::size_t expected_N, std::size_t expected_Ndim,
                           RcpSavedState& s) {
    std::FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) {
        std::cerr << "[rcp_load_state] fopen failed: " << path << "\n";
        return false;
    }
    auto R = [&](void* p, std::size_t n) {
        return std::fread(p, 1, n, fp) == n;
    };
    char magic[4];
    std::uint32_t version = 0;
    std::uint64_t N64 = 0, Ndim64 = 0;
    bool ok = true;
    ok &= R(magic, 4);
    if (std::memcmp(magic, "RCPS", 4) != 0) {
        std::cerr << "[rcp_load_state] bad magic\n";
        std::fclose(fp);
        return false;
    }
    ok &= R(&version, 4);
    if (version != 2) {
        std::cerr << "[rcp_load_state] unsupported version " << version
                  << " (expected 2)\n";
        std::fclose(fp);
        return false;
    }
    ok &= R(&N64, 8);
    ok &= R(&Ndim64, 8);
    ok &= R(&s.step, 8);
    ok &= R(&s.hist_used, 8);
    if (expected_N != 0 && N64 != expected_N) {
        std::cerr << "[rcp_load_state] N mismatch (saved " << N64
                  << " vs expected " << expected_N << ")\n";
        std::fclose(fp);
        return false;
    }
    if (expected_Ndim != 0 && Ndim64 != expected_Ndim) {
        std::cerr << "[rcp_load_state] Ndim mismatch (saved " << Ndim64
                  << " vs expected " << expected_Ndim << ")\n";
        std::fclose(fp);
        return false;
    }
    std::size_t N = N64, Ndim = Ndim64;
    ok &= R(&s.alpha, 8);
    ok &= R(&s.mu, 8);
    ok &= R(&s.kappa, 8);
    ok &= R(&s.m_kappa, 8);
    ok &= R(&s.v_kappa, 8);
    ok &= R(&s.v_update_kappa, 8);
    ok &= R(&s.m_hat_kappa, 8);
    ok &= R(&s.v_hat_kappa, 8);
    ok &= R(&s.t_adam, 8);
    ok &= R(&s.phi, 8);
    ok &= R(&s.phi_modifier, 8);
    ok &= R(&s.phi_max_value, 8);
    ok &= R(&s.phi_max_step, 8);
    ok &= R(&s.mu_change, 8);
    ok &= R(&s.LastPhiUpdate, 8);
    ok &= R(&s.LastAlphaUpdate, 8);
    std::int32_t mu_flag32 = 0, pad32 = 0;
    ok &= R(&mu_flag32, 4);
    ok &= R(&pad32, 4);
    s.mu_flag = mu_flag32;
    ok &= R(s.box, 8 * 8);
    s.positions.assign(N * Ndim, 0.0);
    s.diameters.assign(N, 0.0);
    s.D0.assign(N, 0.0);
    s.m_flat.assign(N * Ndim, 0.0);
    s.v_flat.assign(N * Ndim, 0.0);
    s.phi_history.assign(s.hist_used, 0.0);
    s.F_history.assign(s.hist_used, 0.0);
    s.U_history.assign(s.hist_used, 0.0);
    auto Rvec = [&](std::vector<double>& v) {
        return v.empty() ? true : R(v.data(), v.size() * sizeof(double));
    };
    ok &= Rvec(s.positions);
    ok &= Rvec(s.diameters);
    ok &= Rvec(s.D0);
    ok &= Rvec(s.m_flat);
    ok &= Rvec(s.v_flat);
    ok &= Rvec(s.phi_history);
    ok &= Rvec(s.F_history);
    ok &= Rvec(s.U_history);
    std::fclose(fp);
    if (!ok) {
        std::cerr << "[rcp_load_state] read failed: " << path << "\n";
    } else {
        std::cerr << "[rcp_load_state] loaded step=" << s.step
                  << " hist_used=" << s.hist_used
                  << " mu_flag=" << (int)s.mu_flag
                  << " mu=" << s.mu
                  << " alpha=" << s.alpha << "\n";
    }
    return ok;
}

std::size_t compute_num_changes(const std::vector<std::uint32_t>& refresh) {
    std::size_t sum = 0;
    for (auto v : refresh)
        sum += v;
    return sum;
}

std::vector<std::size_t> sort_indices_by_column(
    const std::vector<std::vector<double>>& x,
    std::size_t col)
{
    std::size_t N = x.size();
    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    auto less = [&](std::size_t a, std::size_t b) {
        double xa = x[a][col];
        double xb = x[b][col];
        if (xa < xb) return true;
        if (xa > xb) return false;
        return a < b;
    };
#if RCP_HAVE_GNU_PARALLEL_SORT
    __gnu_parallel::sort(idx.begin(), idx.end(), less);
#else
    std::sort(idx.begin(), idx.end(), less);
#endif
    return idx;
}

// Cycle 14: same logic but on a flat x[N*Ndim] storage. col indexes the
// dimension to sort by.
std::vector<std::size_t> sort_indices_by_column_flat(
    const double* x,
    std::size_t N,
    std::size_t Ndim,
    std::size_t col)
{
    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    auto less = [&](std::size_t a, std::size_t b) {
        double xa = x[a * Ndim + col];
        double xb = x[b * Ndim + col];
        if (xa < xb) return true;
        if (xa > xb) return false;
        return a < b;
    };
#if RCP_HAVE_GNU_PARALLEL_SORT
    __gnu_parallel::sort(idx.begin(), idx.end(), less);
#else
    std::sort(idx.begin(), idx.end(), less);
#endif
    return idx;
}

// Cycle 8 H17: 3D Morton code (Z-order curve) for spatial locality reorder.
// Each axis quantized to 21 bits → 63-bit Morton code. Particles close in
// 3D space are close in the linear ordering, giving the leaf scan and
// force loop better cache locality on x[]/D[] reads.
static inline std::uint64_t rcp_morton_3d_split21(std::uint64_t v) {
    v &= 0x1FFFFFull;
    v = (v | (v << 32)) & 0x001F00000000FFFFull;
    v = (v | (v << 16)) & 0x001F0000FF0000FFull;
    v = (v | (v <<  8)) & 0x100F00F00F00F00Full;
    v = (v | (v <<  4)) & 0x10C30C30C30C30C3ull;
    v = (v | (v <<  2)) & 0x1249249249249249ull;
    return v;
}
static inline std::uint64_t rcp_morton_3d_encode(std::uint32_t x,
                                                  std::uint32_t y,
                                                  std::uint32_t z) {
    return rcp_morton_3d_split21(x)
         | (rcp_morton_3d_split21(y) << 1)
         | (rcp_morton_3d_split21(z) << 2);
}
std::vector<std::size_t> sort_indices_by_morton_3d(
    const double* x,
    std::size_t N,
    std::size_t Ndim,
    const double* box)
{
    if (Ndim != 3) {
        return sort_indices_by_column_flat(x, N, Ndim, 0);
    }
    const std::uint32_t bits = 21;
    const double scale_x = (double(1u << bits) - 1.0) / box[0];
    const double scale_y = (double(1u << bits) - 1.0) / box[1];
    const double scale_z = (double(1u << bits) - 1.0) / box[2];
    std::vector<std::pair<std::uint64_t, std::uint32_t>> codes(N);
    for (std::size_t i = 0; i < N; ++i) {
        // Clamp positions to [0, box) in case of FP edge cases at boundary.
        double xi = std::max(0.0, std::min(box[0] - 1e-12, x[i*3+0]));
        double yi = std::max(0.0, std::min(box[1] - 1e-12, x[i*3+1]));
        double zi = std::max(0.0, std::min(box[2] - 1e-12, x[i*3+2]));
        std::uint32_t ix = std::uint32_t(xi * scale_x);
        std::uint32_t iy = std::uint32_t(yi * scale_y);
        std::uint32_t iz = std::uint32_t(zi * scale_z);
        codes[i] = {rcp_morton_3d_encode(ix, iy, iz),
                    static_cast<std::uint32_t>(i)};
    }
#if RCP_HAVE_GNU_PARALLEL_SORT
    __gnu_parallel::sort(codes.begin(), codes.end());
#else
    std::sort(codes.begin(), codes.end());
#endif
    std::vector<std::size_t> idx(N);
    for (std::size_t i = 0; i < N; ++i) idx[i] = codes[i].second;
    return idx;
}

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim,
    bool fix_height)
{
    double V_box = 1.0;
    for (double L : box) {
        V_box *= L;
    }

    double sumVol = 0.0;
    for (double Di : D) {
        double radius = Di * 0.5;
        sumVol += sphere_volume(radius, Ndim);
    }

    double factor = 1.0;
    if (fix_height)
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim - 1));
    }
    else
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim));
    }

    std::vector<double> D_scaled = D;
    for (double& Di : D_scaled) {
        Di *= factor;
    }

    return {D_scaled, factor};
}

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
    std::uint32_t* out_max_observed,
    std::uint32_t* out_min_observed)
{
    double Dmax = *std::max_element(D.begin(), D.end());
    double Dmin = *std::min_element(D.begin(), D.end());
    // Median via nth_element — O(N) avg vs O(N log N) for full sort. The
    // downstream `t` calculation depends only on Dsorted[N/2], not the rest
    // of the order, so the partial-sort is equivalent.
    std::vector<double> Dsorted = D;
    std::nth_element(Dsorted.begin(), Dsorted.begin() + N / 2, Dsorted.end());
    double median = Dsorted[N / 2];
    double t = std::min(std::max(median * 1.02, Dmin * 2.25), Dmax) * 0.5;

    const bool probe_sort = rcp_cycle_probe_enabled();
    std::uint64_t sort_tsc0 = probe_sort ? __rdtsc() : 0;
    auto sort_idx = sort_indices_by_column_flat(x, N, Ndim, 0);
    if (probe_sort) {
        g_gpairs_sort_cycles.fetch_add(__rdtsc() - sort_tsc0,
                                       std::memory_order_relaxed);
        g_gpairs_sort_calls.fetch_add(1, std::memory_order_relaxed);
    }
    std::vector<std::size_t> sort_loc(N);
    for (std::size_t k = 0; k < N; ++k)
        sort_loc[sort_idx[k]] = k;

    // Cycle 15d: precompute 1/box[d] so the per-pair MIC wrap uses a
    // multiply instead of a hardware FP division (which is ~14 cyc, not
    // pipelined). Same pattern that Cycle 9 introduced in forces.loop;
    // never propagated to get_pairs until now.
    double gp_inv_box[8] = {0.0};
    for (std::size_t d = 0; d < Ndim; ++d) gp_inv_box[d] = 1.0 / box[d];

    // Cycle 25: build (j, x[j*Ndim+0]) packed array for cell-walk.
    // Measured per-x-load cost in this code was ~2.4 cyc per cand; cell-walk
    // fail-path (90% of cands) only needs sort_idx[pos] and x[j*Ndim+0],
    // so packing them eliminates one cache line load per fail-cand.
    static std::vector<PackedSortEntry> packed_sort_storage;
    if (packed_sort_storage.size() < N) packed_sort_storage.resize(N);
    #pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t jj = sort_idx[k];
        packed_sort_storage[k].j        = static_cast<std::uint32_t>(jj);
        packed_sort_storage[k].D_half_f = static_cast<float>(D[jj] * 0.5);
        packed_sort_storage[k].x0       = x[jj * Ndim];
    }
    const PackedSortEntry* const packed_sort_data = packed_sort_storage.data();

    // Per-particle locks for Phase B cross-writes.
    static std::vector<omp_lock_t> j_locks;
    static std::size_t j_locks_N = 0;
    if (j_locks_N != N) {
        for (std::size_t k = 0; k < j_locks_N; ++k) omp_destroy_lock(&j_locks[k]);
        j_locks.resize(N);
        for (std::size_t k = 0; k < N; ++k) omp_init_lock(&j_locks[k]);
        j_locks_N = N;
    }

    bool warningIssued = false;

    // Cycle 13: flat pairs storage. pairs[i] is now at pairs_data + pair_offsets[i].
    // pairs_data[pair_offsets[i]] is the count; following entries are neighbors.
    // Per-particle cap is max_neighbors_per_particle[i].
    //
    // Algorithm carries Cycle 6's asymmetric scheme: for each refreshed i, walk
    // outward in x-sorted order. For each neighbor j:
    //   - If j > i: append to particle i's slot (lock-free, own data).
    //   - If j < i and refresh[j] == 0: cross-write i into particle j's slot
    //     (per-j lock + dup-check + insert).
    //   - If j < i and refresh[j] == 1: j's own walk handles it.
    const bool probe_gp = rcp_cycle_probe_enabled();
    #pragma omp parallel shared(warningIssued, probe_gp, packed_sort_data, \
                                g_gpairs_cycles, g_gpairs_candidates, g_gpairs_accepted, \
                                g_gpairs_cell_cycles, g_gpairs_dist_cycles, g_gpairs_write_cycles, \
                                g_gpairs_thread_cycles)
    {
        const int gp_tid = omp_get_thread_num();
        std::uint64_t tls_cycles = 0;
        std::uint64_t tls_candidates = 0;
        std::uint64_t tls_accepted = 0;
        // Cycle 15c: split per-candidate cycles into cell-walk (always),
        // dist-test (when first-dim passes), accept-write (when overall passes).
        std::uint64_t tls_cyc_cell = 0;
        std::uint64_t tls_cyc_dist = 0;
        std::uint64_t tls_cyc_write = 0;
        std::uint64_t tsc_start = probe_gp ? __rdtsc() : 0;
        #pragma omp for schedule(dynamic, 32)
        for (std::size_t i = 0; i < N; ++i) {
        if (refresh[i] != 1) continue;
        std::uint32_t* slot_i = pairs_data + pair_offsets[i];
        slot_i[0] = 0;
        const std::uint32_t cap_i = max_neighbors_per_particle[i];

        double r_c_max = (D[i] + 1.01 * Dmax) / 2 + t;
        // Cycle 15j: hoist Di-invariants out of the inner pair loop.
        // Replaces per-candidate `(D[i] + D[j]) / 2 + t` (1 add + 1 div + 1 add)
        // with `Di_half_plus_t + D[j]*0.5` (1 mul + 1 add). Saves ~2 cyc/cand.
        const double Di_half_plus_t = D[i] * 0.5 + t;
        // Cycle 16a: in dense regime, get_pairs is 96% of wall and the
        // inner while loop runs ~700 candidates per particle. Loop-invariants
        // in i go here once; compiler should have done this but being
        // explicit removes any doubt and frees up registers in the hot path.
        const std::size_t ibase = i * Ndim;
        const double xi0 = x[ibase];
        const double xi1 = (Ndim >= 2) ? x[ibase + 1] : 0.0;
        const double xi2 = (Ndim >= 3) ? x[ibase + 2] : 0.0;
        for (int dir : {-1, 1}) {
            std::size_t jdx = 0;
            std::ptrdiff_t delta = 0;
            bool go = true;
            while (go) {
                ++jdx;
                delta += dir;
                if (jdx > N / 2) { go = false; break; }

                // Cycle 15g: replace `% N` (hardware integer divide, ~25 cyc,
                // not pipelined) with two predicted branches. Invariant
                // jdx <= N/2 means |dir*jdx| < N, so a single +/- N suffices
                // to wrap. Same modular arithmetic, ~10x cheaper per cand.
                // Cycle 18: incremental signed_pos. `dir * jdx` recomputed
                // each iter is now `delta += dir` carried by the loop. Saves
                // the multiply per cand and lets the compiler keep the
                // wrap as a single conditional sub/add.
                std::ptrdiff_t signed_pos =
                    static_cast<std::ptrdiff_t>(sort_loc[i]) + delta;
                if (signed_pos < 0)
                    signed_pos += static_cast<std::ptrdiff_t>(N);
                else if (signed_pos >= static_cast<std::ptrdiff_t>(N))
                    signed_pos -= static_cast<std::ptrdiff_t>(N);
                std::size_t pos = static_cast<std::size_t>(signed_pos);
                // Cycle 25: read packed (j, x[j*Ndim+0]) in a single 16-byte
                // load. Eliminates the separate x[jbase] cache line load on
                // the 90% fail-path. The pass-path (~10%) still loads
                // x[jbase+1, +2] for d2 below — different cache line now,
                // small net cost — but the average wins by far.
                const PackedSortEntry& pe = packed_sort_data[pos];
                const std::size_t j = pe.j;
                double dx = pe.x0 - xi0;
                dx -= std::nearbyint(dx * gp_inv_box[0]) * box[0];

                bool first_dim_pass = (std::fabs(dx) <= r_c_max);
                if (!first_dim_pass) { go = false; break; }

                // Cycle 26: r_c uses pre-packed D[j]*0.5 from the same
                // cache line as j and x0. Avoids the separate D[j] load.
                double r_c = Di_half_plus_t + static_cast<double>(pe.D_half_f);
                if (std::fabs(dx) < r_c) {
                    // Cycle 25: jbase computed here (only needed on pass-path).
                    const std::size_t jbase = j * Ndim;
                    double d2 = dx * dx;
                    if (Ndim == 3) {
                        double dy = x[jbase + 1] - xi1;
                        dy -= std::nearbyint(dy * gp_inv_box[1]) * box[1];
                        double dz = x[jbase + 2] - xi2;
                        dz -= std::nearbyint(dz * gp_inv_box[2]) * box[2];
                        d2 += dy * dy + dz * dz;
                    } else {
                        for (std::size_t d = 1; d < Ndim; ++d) {
                            double dz = x[jbase + d] - x[ibase + d];
                            dz -= std::nearbyint(dz * gp_inv_box[d]) * box[d];
                            d2 += dz * dz;
                        }
                    }
                    if (d2 >= r_c * r_c) continue;

                    if (j > i) {
                        // Lockless append to particle i's slot.
                        std::uint32_t& cnt = slot_i[0];
                        if (cnt < cap_i) {
                            ++cnt;
                            slot_i[cnt] = static_cast<std::uint32_t>(j);
                        } else {
                            #pragma omp critical(warn_neighbor_overflow_A)
                            if (!warningIssued) {
                                std::cerr << "Error: particle " << i
                                          << " exceeded MaxNb_i = " << cap_i
                                          << " (size-based per-particle cap). "
                                          << "Raise MaxNbBig in the run options.\n";
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    } else if (refresh[j] == 0) {
                        // Cross-write to particle j's slot.
                        std::uint32_t* slot_j = pairs_data + pair_offsets[j];
                        const std::uint32_t cap_j = max_neighbors_per_particle[j];
                        omp_set_lock(&j_locks[j]);
                        bool duplicate = false;
                        std::uint32_t cntj = slot_j[0];
                        for (std::uint32_t kk = 1; kk <= cntj; ++kk) {
                            if (slot_j[kk] == i) { duplicate = true; break; }
                        }
                        if (!duplicate) {
                            if (cntj < cap_j) {
                                slot_j[0] = cntj + 1;
                                slot_j[cntj + 1] = static_cast<std::uint32_t>(i);
                                omp_unset_lock(&j_locks[j]);
                            } else {
                                omp_unset_lock(&j_locks[j]);
                                #pragma omp critical(warn_neighbor_overflow_B)
                                if (!warningIssued) {
                                    std::cerr << "Error: particle " << j
                                              << " exceeded MaxNb_j = " << cap_j
                                              << " (size-based per-particle cap). "
                                              << "Raise MaxNbBig in the run options.\n";
                                    std::exit(EXIT_FAILURE);
                                }
                            }
                        } else {
                            omp_unset_lock(&j_locks[j]);
                        }
                    }
                }
            }
        }
        }  // close omp for body

        if (probe_gp) {
            tls_cycles = __rdtsc() - tsc_start;
            g_gpairs_cycles.fetch_add(tls_cycles, std::memory_order_relaxed);
            g_gpairs_candidates.fetch_add(tls_candidates, std::memory_order_relaxed);
            g_gpairs_accepted.fetch_add(tls_accepted, std::memory_order_relaxed);
            g_gpairs_cell_cycles.fetch_add(tls_cyc_cell, std::memory_order_relaxed);
            g_gpairs_dist_cycles.fetch_add(tls_cyc_dist, std::memory_order_relaxed);
            g_gpairs_write_cycles.fetch_add(tls_cyc_write, std::memory_order_relaxed);
            if (gp_tid >= 0 && gp_tid < RCP_PROBE_MAX_THREADS) {
                g_gpairs_thread_cycles[gp_tid].fetch_add(tls_cycles, std::memory_order_relaxed);
            }
        }
    }  // close omp parallel block

    // Canonicalize: sort each particle's neighbor list.
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        std::uint32_t* slot = pairs_data + pair_offsets[i];
        std::uint32_t cnt = slot[0];
        if (cnt > 1) {
            std::sort(slot + 1, slot + 1 + cnt);
        }
    }

    // Cycle 13 diagnostics: track max/min observed neighbor count.
    if (out_max_observed || out_min_observed) {
        std::uint32_t local_max = 0;
        std::uint32_t local_min = std::numeric_limits<std::uint32_t>::max();
        #pragma omp parallel for schedule(static) \
            reduction(max:local_max) reduction(min:local_min)
        for (std::size_t i = 0; i < N; ++i) {
            std::uint32_t cnt = pairs_data[pair_offsets[i]];
            if (cnt > local_max) local_max = cnt;
            if (cnt < local_min) local_min = cnt;
        }
        if (out_max_observed) *out_max_observed = local_max;
        if (out_min_observed) *out_min_observed = local_min;
    }

    ++g_pairs_version;
    (void)walls;
}

// ============================================================================
// KDTREE PROJECT — Cycle 1 scaffolding (namespace-scope: kdtree backend +
// dispatch). File-scope selector and shadow counters live earlier in the file.
// ============================================================================

// Per-particle set-equality check between two pairs_data buffers using the
// same pair_offsets layout. Updates global counters.
static void compare_pairs_sets(
    std::size_t N,
    const std::uint64_t* pair_offsets,
    const std::uint32_t* pairs_oracle,
    const std::uint32_t* pairs_test)
{
    // Optional verbose mode: set RCP_SHADOW_VERBOSE=1 to print the first few
    // divergent particles (with missing / extra j lists) at the first
    // divergent steps. Helpful for kd-tree bring-up; harmless when off.
    static int verbose = []() {
        const char* s = std::getenv("RCP_SHADOW_VERBOSE");
        return (s && *s && *s != '0') ? 1 : 0;
    }();
    std::uint64_t local_particles_diff = 0;
    std::uint64_t local_missing = 0;
    std::uint64_t local_extra = 0;
    std::uint64_t first_local = UINT64_MAX;
    int printed = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const std::uint32_t* slot_o = pairs_oracle + pair_offsets[i];
        const std::uint32_t* slot_t = pairs_test   + pair_offsets[i];
        const std::uint32_t cnt_o = slot_o[0];
        const std::uint32_t cnt_t = slot_t[0];
        std::set<std::uint32_t> set_o(slot_o + 1, slot_o + 1 + cnt_o);
        std::set<std::uint32_t> set_t(slot_t + 1, slot_t + 1 + cnt_t);
        if (set_o != set_t) {
            ++local_particles_diff;
            for (auto j : set_o) if (!set_t.count(j)) ++local_missing;
            for (auto j : set_t) if (!set_o.count(j)) ++local_extra;
            if (first_local == UINT64_MAX) first_local = static_cast<std::uint64_t>(i);
            if (verbose && printed < 3) {
                std::cerr << "[shadow diff] step=" << g_shadow_current_step
                          << " particle=" << i
                          << " oracle_cnt=" << cnt_o << " test_cnt=" << cnt_t << "\n";
                std::cerr << "    missing in test:";
                for (auto j : set_o) if (!set_t.count(j)) std::cerr << " " << j;
                std::cerr << "\n    extra in test:  ";
                for (auto j : set_t) if (!set_o.count(j)) std::cerr << " " << j;
                std::cerr << "\n";
                ++printed;
            }
        }
    }
    g_shadow_calls.fetch_add(1, std::memory_order_relaxed);
    if (local_particles_diff) {
        g_shadow_particles_diff.fetch_add(local_particles_diff, std::memory_order_relaxed);
        g_shadow_pairs_missing_in_kdtree.fetch_add(local_missing, std::memory_order_relaxed);
        g_shadow_pairs_extra_in_kdtree.fetch_add(local_extra, std::memory_order_relaxed);
        std::uint64_t expected = UINT64_MAX;
        if (g_shadow_first_diff_step.compare_exchange_strong(
                expected, g_shadow_current_step, std::memory_order_relaxed)) {
            g_shadow_first_diff_particle.store(first_local, std::memory_order_relaxed);
        }
    }
}

// ============================================================================
// Cycle 2: real kd-tree backend implementation
// ============================================================================
// Design: simple recursive kd-tree built fresh per call. Median split along
// widest dimension. Leaf size = 8. Range queries enumerate PBC images per
// axis (typically 1, sometimes 2 when near box boundary).
//
// Correctness invariants (mirror get_pairs_nd_3 exactly):
//   - same shell t formula
//   - same r_c_max(i) = (D_i + 1.01·Dmax)/2 + t
//   - same inclusion: dist_PBC(i,j) < (D_i + D_j)/2 + t
//   - same asymmetric storage: j > i own-write, j < i cross-write with refresh-gated lock
//   - same per-particle MaxNb_i cap with same error on overflow
//   - same shared j_locks (declared static inside get_pairs_nd_3; we reuse a
//     copy here — Cycle 2 lives or dies on logic match, not lock-sharing)

namespace {

// kdtree Cycle 4: KD_MAX_NDIM reduced from 8 -> 4. KdNode bbox arrays shrink
// from 64 -> 32 bytes each, total per-node footprint 144 -> 80 bytes. For
// N=100K that takes the tree from ~7.2 MB to ~4 MB — fits L3 comfortably.
// Production runs are Ndim=2 or 3; supporting Ndim=4 covers the realistic
// extensions. For Ndim>4 the user must increase this constant and rebuild.
constexpr std::size_t KD_MAX_NDIM = 4;
// kdtree Cycle 4: LEAF_SIZE 8 -> 64. Larger leaves shorten traversal depth
// (three fewer levels for N=100K), shrink the node count, and trade
// tree-recursion overhead for in-leaf scan work — the in-leaf scan is fully
// unrolled for Ndim=3 so each extra particle is just a few cycles. Swept
// LEAF in {8,16,32,64,128} at H6: 64 was the empirical sweet spot.
// Re-tested at H11 (post-deferred-cross-write): LEAF=128 → 237.58 s
// (regression vs 64=218.83 s) because larger leaves visit more particles
// per scan (most rejected by per-pair cutoff), outweighing reduced
// traversal depth. 64 remains optimal.
constexpr std::uint32_t KD_LEAF_SIZE = 32;
constexpr std::uint32_t KD_NULL = 0xFFFFFFFFu;

// Cycle 8 H18: KdNode compacted to 44 bytes (down from ~88) so the tree
// fits L2 even at N=500K (1.4MB vs 2.5MB). float bbox stored with proper
// rounding direction (lo→down, hi→up, max_D→up) so the float bbox is a
// SUPERSET of the true double bbox — pruning stays conservative, no
// missed neighbors. The outer MIC test in get_pairs_kdtree's accept-write
// path runs on double positions, so set-equivalence with sort+walk holds.
// Restricts the kdtree backend to Ndim ≤ 3 (production constraint).
struct KdNode {
    std::uint32_t children[2] = {KD_NULL, KD_NULL};  // 8 B
    std::uint32_t first = 0;                          // 4 B
    std::uint32_t count = 0;                          // 4 B
    float bbox_lo[3]{};                               // 12 B
    float bbox_hi[3]{};                               // 12 B
    float max_D = 0.0f;                               // 4 B
    // Total: 44 B (single cache line holds ~1.45 nodes)
};
static_assert(sizeof(KdNode) == 44,
              "KdNode should be 44 bytes after H18 compaction");

// Cycle 10 H43: octree node. 8-way branching cuts tree depth to log_8(N/LEAF)
// (~5 for N=500K LEAF=32 vs binary's 14). Each internal node has 8 children
// corresponding to spatial octants. Per-query nodes visited drops from
// ~30 (binary) to ~12 (octree), and the 8 children's bbox checks can be
// SIMD-batched (AVX2 4-wide).
//
// Size: 32 (children) + 4 (first) + 4 (count) + 12 (bbox_lo) + 12 (bbox_hi)
// + 4 (max_D) = 68 bytes. Larger per node but fewer nodes overall (octree
// has ~N/LEAF*8/7 nodes vs binary's ~2*N/LEAF, so total memory is comparable).
struct OctNode {
    std::uint32_t children[8] = {KD_NULL, KD_NULL, KD_NULL, KD_NULL,
                                  KD_NULL, KD_NULL, KD_NULL, KD_NULL};
    std::uint32_t first = 0;
    std::uint32_t count = 0;
    float bbox_lo[3]{};
    float bbox_hi[3]{};
    float max_D = 0.0f;
};
static_assert(sizeof(OctNode) == 68,
              "OctNode should be 68 bytes after H43");

// Cycle 8 H18: directional float rounding so bbox/max_D rounded to float
// remains a conservative bound on the true double values. No missed
// neighbors during tree pruning.
static inline float rcp_f_round_down(double v) {
    float f = static_cast<float>(v);
    if (static_cast<double>(f) > v)
        f = std::nextafterf(f, -std::numeric_limits<float>::infinity());
    return f;
}
static inline float rcp_f_round_up(double v) {
    float f = static_cast<float>(v);
    if (static_cast<double>(f) < v)
        f = std::nextafterf(f, std::numeric_limits<float>::infinity());
    return f;
}

struct KdTree {
    std::vector<KdNode> nodes;
    std::vector<std::uint32_t> perm;  // particle indices in tree-leaf order
    std::atomic<std::uint32_t> next_node{0};  // kdtree Cycle 4: parallel-build allocator
    bool g_kd_root_pre_partitioned = false;  // Cycle 8 H19
    std::uint32_t g_kd_root_mid = 0;  // root split point after partition
    // kdtree Cycle 6 H2: incremental-patch infrastructure. parent_of_node lets
    // us walk from a leaf to root and expand ancestor bboxes when a particle
    // moves; leaf_of_particle maps each particle index to the leaf node it
    // lives in so patching is O(log N) per moved particle.
    std::vector<std::uint32_t> parent_of_node;
    std::vector<std::uint32_t> leaf_of_particle;
    // Cycle 10 H30: per-node depth used by the dedup'd bottom-up propagate
    // path. Set in build_recursive; uint8 is plenty (max depth ≈ 32 for any
    // realistic N).
    std::vector<std::uint8_t> depth_of_node;

    // Build from particle positions x (length N*Ndim). All particles assumed
    // already inside the canonical box [0, box[d]) for each d.
    //
    // kdtree Cycle 4: parallel recursive build.
    // kdtree Cycle 6 H3: D_arr is per-particle diameter; used to populate
    // node.max_D bottom-up so the range-query traversal can prune subtrees
    // that contain no large particles tighter than the global (Di+Dmax)/2+t.
    void build(const double* x, const double* D_arr,
               std::size_t N, std::size_t Ndim) {
        if (N == 0) {
            nodes.clear();
            perm.clear();
            next_node.store(0, std::memory_order_relaxed);
            return;
        }
        // perm: fill 0..N-1 (overwrites every entry, no need for assign+zero)
        if (perm.size() < N) perm.resize(N);
        for (std::size_t k = 0; k < N; ++k)
            perm[k] = static_cast<std::uint32_t>(k);
        // Upper bound on node count. Median splits keep the tree roughly
        // balanced: a parent of count > LEAF splits into two halves of
        // count/2 each, recursion stops at count <= LEAF, so leaves end up
        // with at least floor((LEAF+1)/2) and at most LEAF particles. Number
        // of leaves <= 2*N/LEAF and total nodes <= 4*N/LEAF. Plus a small
        // padding for the root-level edge cases. Keeps the tree small enough
        // to fit in L3 for production sizes; pre-sized so atomic-counter
        // writes from concurrent build tasks never race against a realloc.
        // build_recursive writes every field of every reachable node, so no
        // zero-init is required — just resize-grow on first build (or when
        // N grows) and reuse the buffer thereafter.
        const std::size_t upper = 4 * N / KD_LEAF_SIZE + 64;
        if (nodes.size() < upper) nodes.resize(upper);
        if (parent_of_node.size() < upper) parent_of_node.resize(upper, KD_NULL);
        if (leaf_of_particle.size() < N) leaf_of_particle.resize(N);
        if (depth_of_node.size() < upper) depth_of_node.resize(upper, 0);
        // Cycle 8 H20 (packed-position leaf storage) tried + reverted.
        parent_of_node[0] = KD_NULL;  // root has no parent
        next_node.store(1, std::memory_order_relaxed);  // root at index 0
        // Cycle 8 H19: parallel root partition via __gnu_parallel.
        // (Cycle 9 H22 tested manual sample-based parallel partition with
        // local less/greater buffers; performance was within noise of H19,
        // so reverted to simpler __gnu_parallel approach.)
#if RCP_HAVE_GNU_PARALLEL_SORT
        if (N > KD_LEAF_SIZE) {
            const std::size_t split_dim_root = 0;
            const std::uint32_t mid = static_cast<std::uint32_t>(N) / 2;
            __gnu_parallel::nth_element(
                perm.begin(), perm.begin() + mid, perm.end(),
                [x, Ndim, split_dim_root](std::uint32_t a, std::uint32_t b) {
                    return x[a * Ndim + split_dim_root]
                         < x[b * Ndim + split_dim_root];
                });
            g_kd_root_mid = mid;
            g_kd_root_pre_partitioned = true;
        }
#endif
        #pragma omp parallel
        {
            #pragma omp single
            build_recursive(x, D_arr, Ndim, 0, 0,
                            static_cast<std::uint32_t>(N), 0);
            // implicit barrier — all spawned tasks complete here
        }
        g_kd_root_pre_partitioned = false;
        // Cycle 10 H34 (BFS reorder of node array) tested + reverted —
        // wall save was 0 ± noise at N=500K. Cache-miss model overestimated
        // tree-traversal cost; tree (~1.4MB nodes) fits L2 even in DFS
        // layout, so reorder provided no measurable benefit while adding
        // ~5s of per-build overhead.
    }

    // kdtree Cycle 6 H2v2: TRUE bbox-recompute incremental patch. Unlike
    // expand-only patch (H2v1), this recomputes the leaf's bbox tightly
    // from current positions of all its particles, then walks up to root
    // recomputing each ancestor from its two children's bboxes. Result:
    // tree's bboxes stay accurate (no bloat) for the modified leaves'
    // ancestor chain. Unmodified leaves keep stale bboxes — covered by
    // a Dmin/2 query-radius inflation handled by the caller.
    void recompute_leaf_bbox(std::uint32_t leaf_idx,
                             const double* x, const double* D_arr,
                             std::size_t Ndim) {
        KdNode& L = nodes[leaf_idx];
        // Cycle 8 H18: compute in double, store as conservatively-rounded
        // float so the float bbox is a SUPERSET of the true double bbox.
        double lo[3] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[3] = {-std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
        double mx_D = 0.0;
        const std::uint32_t end = L.first + L.count;
        for (std::uint32_t k = L.first; k < end; ++k) {
            const std::uint32_t p = perm[k];
            const double* xp = x + p * Ndim;
            for (std::size_t d = 0; d < Ndim; ++d) {
                if (xp[d] < lo[d]) lo[d] = xp[d];
                if (xp[d] > hi[d]) hi[d] = xp[d];
            }
            const double Dk = D_arr[p];
            if (Dk > mx_D) mx_D = Dk;
        }
        for (std::size_t d = 0; d < Ndim; ++d) {
            L.bbox_lo[d] = rcp_f_round_down(lo[d]);
            L.bbox_hi[d] = rcp_f_round_up(hi[d]);
        }
        L.max_D = rcp_f_round_up(mx_D);
    }

    // Walk leaf→root, recomputing each ancestor's bbox/max_D as the
    // enclosing of its two children's CURRENT values. Children's float
    // bboxes are already conservatively rounded; union remains
    // conservative (min of two lo bounds stays a lo bound, etc.).
    void propagate_up_from(std::uint32_t leaf_idx, std::size_t Ndim) {
        std::uint32_t curr = parent_of_node[leaf_idx];
        while (curr != KD_NULL) {
            KdNode& n = nodes[curr];
            const KdNode& lc = nodes[n.children[0]];
            const KdNode& rc = nodes[n.children[1]];
            for (std::size_t d = 0; d < Ndim; ++d) {
                n.bbox_lo[d] = std::min(lc.bbox_lo[d], rc.bbox_lo[d]);
                n.bbox_hi[d] = std::max(lc.bbox_hi[d], rc.bbox_hi[d]);
            }
            n.max_D = std::max(lc.max_D, rc.max_D);
            curr = parent_of_node[curr];
        }
    }

    // Sphere range query with per-pair inclusion cutoff (Cycle 6 H3).
    // For each particle p, the cutoff is (Di + D[p])/2 + t_shell — we only
    // include p if its Euclidean distance to xq is below this. The
    // traversal also uses each subtree's max_D to tighten bbox pruning to
    // (Di + node.max_D)/2 + t_shell, much smaller than the global
    // (Di + Dmax)/2 + t for subtrees containing only small particles.
    void range_query(const double* xq, double Di, double t_shell,
                     const double* x, const double* D_arr,
                     std::size_t Ndim,
                     std::vector<std::uint32_t>& out) const {
        if (nodes.empty()) return;
        if (Ndim == 3) {
            range_query_recursive_3d(0, xq, Di, t_shell, x, D_arr, out);
        } else {
            range_query_recursive(0, xq, Di, t_shell, x, D_arr, Ndim, out);
        }
    }

private:
    // Below this subtree count, recurse serially (task overhead > work).
    static constexpr std::uint32_t KD_PARALLEL_SPLIT_MIN = 4096;

    // Build the subtree rooted at node_idx covering perm[begin..end).
    // kdtree Cycle 6 H6: round-robin split dim + bbox-from-children. No
    // per-particle bbox compute at internal nodes (only at leaves, where we
    // need them tight). Internal bboxes are folded up from children's
    // bboxes after recursion. Split dimension cycles by depth (depth%Ndim);
    // for particle distributions that fill the box quasi-uniformly this
    // gives a balanced tree without the O(N) widest-dim probe at each level.
    void build_recursive(const double* x, const double* D_arr,
                         std::size_t Ndim,
                         std::uint32_t node_idx,
                         std::uint32_t begin, std::uint32_t end,
                         std::uint32_t depth) {
        KdNode& node = nodes[node_idx];
        // Cycle 10 H30: record depth for the patch-path bottom-up dedup.
        depth_of_node[node_idx] = static_cast<std::uint8_t>(
            depth > 255u ? 255u : depth);
        const std::uint32_t count = end - begin;
        if (count <= KD_LEAF_SIZE) {
            // Leaf: compute tight bbox from particle positions + max_D.
            // Cycle 8 H18: accumulate in double, store as conservatively
            // rounded float so bbox is a SUPERSET of the true bbox.
            double lo[3] = {std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()};
            double hi[3] = {-std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()};
            double mx_D = 0.0;
            for (std::uint32_t k = begin; k < end; ++k) {
                const std::uint32_t p = perm[k];
                const double* xp = x + p * Ndim;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    const double v = xp[d];
                    if (v < lo[d]) lo[d] = v;
                    if (v > hi[d]) hi[d] = v;
                }
                leaf_of_particle[p] = node_idx;
                const double Dk = D_arr[p];
                if (Dk > mx_D) mx_D = Dk;
            }
            for (std::size_t d = 0; d < Ndim; ++d) {
                node.bbox_lo[d] = rcp_f_round_down(lo[d]);
                node.bbox_hi[d] = rcp_f_round_up(hi[d]);
            }
            node.max_D = rcp_f_round_up(mx_D);
            node.first = begin;
            node.count = count;
            node.children[0] = KD_NULL;
            node.children[1] = KD_NULL;
            return;
        }

        // Internal: round-robin split dim (no bbox probe needed).
        const std::size_t split_dim = depth % Ndim;
        // Cycle 9 H22: at root, use the parallel-partition's actual mid
        // (may differ slightly from begin+count/2 due to sample-based
        // pivot). For deeper levels, exact median split as usual.
        std::uint32_t mid;
        if (depth == 0 && g_kd_root_pre_partitioned) {
            mid = g_kd_root_mid;  // from parallel partition
            // No nth_element — partition is done.
        } else {
            mid = begin + count / 2;
            std::nth_element(
                perm.begin() + begin, perm.begin() + mid, perm.begin() + end,
                [x, Ndim, split_dim](std::uint32_t a, std::uint32_t b) {
                    return x[a * Ndim + split_dim] < x[b * Ndim + split_dim];
                });
        }

        // Allocate child node slots from the atomic counter.
        // Cycle 10 H41: one fetch_add(2) instead of two fetch_add(1)s.
        // Halves atomic ops and guarantees sibling pair lands at
        // consecutive indices even under parallel tasks (interleaved
        // tasks no longer scatter siblings).
        const std::uint32_t left_idx = next_node.fetch_add(2, std::memory_order_relaxed);
        const std::uint32_t right_idx = left_idx + 1;

        node.children[0] = left_idx;
        node.children[1] = right_idx;
        node.count = 0;  // mark internal
        parent_of_node[left_idx]  = node_idx;
        parent_of_node[right_idx] = node_idx;

        if (count >= KD_PARALLEL_SPLIT_MIN) {
            #pragma omp task firstprivate(left_idx, begin, mid) shared(x, D_arr)
            build_recursive(x, D_arr, Ndim, left_idx, begin, mid, depth + 1);
            #pragma omp task firstprivate(right_idx, mid, end) shared(x, D_arr)
            build_recursive(x, D_arr, Ndim, right_idx, mid, end, depth + 1);
            #pragma omp taskwait
        } else {
            build_recursive(x, D_arr, Ndim, left_idx, begin, mid, depth + 1);
            build_recursive(x, D_arr, Ndim, right_idx, mid, end, depth + 1);
        }
        // Cycle 6 H6: fold internal node's bbox and max_D from children.
        const KdNode& lc = nodes[left_idx];
        const KdNode& rc = nodes[right_idx];
        for (std::size_t d = 0; d < Ndim; ++d) {
            node.bbox_lo[d] = std::min(lc.bbox_lo[d], rc.bbox_lo[d]);
            node.bbox_hi[d] = std::max(lc.bbox_hi[d], rc.bbox_hi[d]);
        }
        node.max_D = std::max(lc.max_D, rc.max_D);
    }

    // Min squared distance from point xq to node's bbox (generic Ndim).
    static double bbox_min_d2(const double* xq, const double* lo,
                              const double* hi, std::size_t Ndim) {
        double d2 = 0.0;
        for (std::size_t d = 0; d < Ndim; ++d) {
            double diff = 0.0;
            if (xq[d] < lo[d]) diff = lo[d] - xq[d];
            else if (xq[d] > hi[d]) diff = xq[d] - hi[d];
            d2 += diff * diff;
        }
        return d2;
    }

    // kdtree Cycle 4 + Cycle 6 H3: Ndim=3 fast path with per-subtree pruning
    // by node.max_D and per-pair cutoff in the leaf scan.
    // Cycle 7 H14: iterative range query (explicit stack). Eliminates
    // recursion overhead — for our typical ~30-nodes-per-query traversal,
    // saves ~100 cyc/node × 30 = 3000 cyc/query in function call setup.
    // Prefetch + branchless bbox check carry over from H12/H13.
    //
    // Cycle 8 H16: function-local fast-math enables FP reassociation so
    // the compiler can auto-vectorize the leaf-scan inner loop. The MIC
    // test in get_pairs_kdtree's outer accept-write path is NOT marked
    // fast-math — it stays IEEE-deterministic to preserve byte-equivalence
    // with sort+walk on borderline pairs.
    __attribute__((optimize("fast-math","tree-vectorize")))
    void range_query_recursive_3d(std::uint32_t root_idx, const double* xq,
                                  double Di, double t_shell,
                                  const double* x, const double* D_arr,
                                  std::vector<std::uint32_t>& out) const {
        constexpr std::size_t MAX_STACK = 64;
        std::uint32_t stack[MAX_STACK];
        int top = 0;
        stack[top++] = root_idx;
        const double xq0 = xq[0], xq1 = xq[1], xq2 = xq[2];
        while (top > 0) {
            const std::uint32_t node_idx = stack[--top];
            const KdNode& n = nodes[node_idx];
            // Prefetch children early (H13).
            if (n.children[0] != KD_NULL) {
                __builtin_prefetch(&nodes[n.children[0]], 0, 0);
                __builtin_prefetch(&nodes[n.children[1]], 0, 0);
            }
            // Branchless bbox-min-d2 for Ndim=3 (H12).
            const double diff0 = std::max(0.0,
                std::max(n.bbox_lo[0] - xq0, xq0 - n.bbox_hi[0]));
            const double diff1 = std::max(0.0,
                std::max(n.bbox_lo[1] - xq1, xq1 - n.bbox_hi[1]));
            const double diff2 = std::max(0.0,
                std::max(n.bbox_lo[2] - xq2, xq2 - n.bbox_hi[2]));
            const double dbb2 = diff0*diff0 + diff1*diff1 + diff2*diff2;
            const double r_node = (Di + n.max_D) * 0.5 + t_shell;
            if (dbb2 > r_node * r_node) continue;
            if (n.children[0] == KD_NULL) {
                // Cycle 8 H20: tested packed-position SoA leaf storage —
                // regressed wall by ~5s. The packing work at build/patch
                // time exceeded the scan savings; HW prefetcher already
                // covered the indirect x[perm[k]*3] reads well (post-H17
                // Morton, perm is spatially clustered). Reverted to plain
                // indirect-access scan.
                const std::uint32_t cnt = n.count;
                const std::uint32_t first = n.first;
                for (std::uint32_t k = 0; k < cnt; ++k) {
                    const std::uint32_t p = perm[first + k];
                    const double* xp = x + p * 3u;
                    const double dd0 = xp[0] - xq0;
                    const double dd1 = xp[1] - xq1;
                    const double dd2 = xp[2] - xq2;
                    const double d2 = dd0*dd0 + dd1*dd1 + dd2*dd2;
                    const double r_p = (Di + D_arr[p]) * 0.5 + t_shell;
                    if (d2 <= r_p * r_p) out.push_back(p);
                }
                continue;
            }
            // Push children onto stack.
            stack[top++] = n.children[0];
            stack[top++] = n.children[1];
        }
    }

    void range_query_recursive(std::uint32_t node_idx, const double* xq,
                               double Di, double t_shell,
                               const double* x, const double* D_arr,
                               std::size_t Ndim,
                               std::vector<std::uint32_t>& out) const {
        const KdNode& n = nodes[node_idx];
        const double r_node = (Di + double(n.max_D)) * 0.5 + t_shell;
        // Cycle 8 H18: inline bbox_min_d2 here so we can read the float
        // bbox without converting via a temporary double array.
        double dbb2 = 0.0;
        for (std::size_t d = 0; d < Ndim; ++d) {
            const double lo_d = double(n.bbox_lo[d]);
            const double hi_d = double(n.bbox_hi[d]);
            double diff = 0.0;
            if (xq[d] < lo_d) diff = lo_d - xq[d];
            else if (xq[d] > hi_d) diff = xq[d] - hi_d;
            dbb2 += diff * diff;
        }
        if (dbb2 > r_node * r_node) return;
        if (n.children[0] == KD_NULL) {
            // Leaf
            for (std::uint32_t k = 0; k < n.count; ++k) {
                const std::uint32_t p = perm[n.first + k];
                double d2 = 0.0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    const double dd = x[p * Ndim + d] - xq[d];
                    d2 += dd * dd;
                }
                const double r_p = (Di + D_arr[p]) * 0.5 + t_shell;
                if (d2 <= r_p * r_p) out.push_back(p);
            }
            return;
        }
        range_query_recursive(n.children[0], xq, Di, t_shell, x, D_arr, Ndim, out);
        range_query_recursive(n.children[1], xq, Di, t_shell, x, D_arr, Ndim, out);
    }
};

// ===========================================================================
// Cycle 10 H43: OctTree (8-way branching) — alternative to binary KdTree.
//
// Rationale: range_query at N=500K is compute-bound (~90 cyc/node, mostly
// per-node bbox check), and binary tree visits ~30 nodes/query. Octree
// visits ~12 nodes/query but each internal node check involves 8 children's
// bboxes. Without SIMD, octree compute per query is roughly equivalent to
// binary; with AVX2 4-wide SIMD across children, the per-node compute drops
// 4x, giving ~3x net query speedup.
//
// Phase A (this commit): correctness — scalar serial build, scalar 8-way
// range_query, byte-test against sort+walk oracle via shadow mode.
// Phase B-E: SIMD bbox check, OMP parallel build, patch path integration.
// ===========================================================================
struct OctTree {
    std::vector<OctNode> nodes;
    std::vector<std::uint32_t> perm;  // particle indices in tree-leaf order
    std::atomic<std::uint32_t> next_node{0};
    std::vector<std::uint32_t> parent_of_node;
    std::vector<std::uint32_t> leaf_of_particle;
    std::vector<std::uint8_t> depth_of_node;
    // Cycle 11 H38: sorted Morton codes (one per entry in perm[]). Filled by
    // the LBVH build path; used by build_recursive_lbvh to identify octant
    // boundaries via Morton bit triplets instead of median-split nth_element.
    std::vector<std::uint64_t> sorted_morton;

    void build(const double* x, const double* D_arr,
               std::size_t N, std::size_t Ndim) {
        if (N == 0 || Ndim != 3) {
            nodes.clear();
            perm.clear();
            next_node.store(0, std::memory_order_relaxed);
            return;
        }
        if (perm.size() < N) perm.resize(N);
        for (std::size_t k = 0; k < N; ++k)
            perm[k] = static_cast<std::uint32_t>(k);
        // Upper bound on node count. Octree allocates 8 child slots per
        // internal node, including empty children. For balanced tree of
        // depth D, total slots = (8^(D+1)-1)/7. D = ceil(log_8(N/LEAF)).
        // For N=500K LEAF=32, D=5 → ~37K slots. Conservative bound that
        // handles all N up to 16M: 16 × N/LEAF + 256.
        const std::size_t upper = 16 * (N / KD_LEAF_SIZE + 1) + 256;
        if (nodes.size() < upper) nodes.resize(upper);
        if (parent_of_node.size() < upper)
            parent_of_node.resize(upper, KD_NULL);
        if (depth_of_node.size() < upper) depth_of_node.resize(upper, 0);
        if (leaf_of_particle.size() < N) leaf_of_particle.resize(N);
        parent_of_node[0] = KD_NULL;
        next_node.store(1, std::memory_order_relaxed);
        // Cycle 11 H38: LBVH (Linear-BVH) build via Morton sort. Replaces the
        // nth_element-based recursive partition with a sort-then-scan
        // approach. The dominant cost is the sort, which is parallelizable
        // via __gnu_parallel::sort — no OMP task scheduler pitfalls. The
        // subsequent build_recursive_lbvh is a serial linear scan, O(N) total.
        build_lbvh(x, D_arr, static_cast<std::uint32_t>(N));
    }

    // Cycle 11 H38: LBVH build. Computes particle bbox, encodes Morton 3D
    // per particle, sorts perm[] by Morton via parallel sort, then recurses
    // via Morton-bit-triplet boundary scanning instead of median split.
    void build_lbvh(const double* x, const double* D_arr, std::uint32_t N) {
        // Cycle 12 H48 (parallel build prelude — bbox + Morton + extract)
        // tested + reverted. The 3 OMP regions per build call added ~30us
        // overhead each (~90us per build), and the parallel loops only saved
        // a few ms of work on each. At our call rate this was net-negative.
        // Reverted to serial.
        // 1. Compute particle bbox in serial (cheap O(N)).
        double box_min[3] = {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity()};
        double box_max[3] = {-std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()};
        for (std::uint32_t i = 0; i < N; ++i) {
            const double* xi = x + i * 3u;
            for (int d = 0; d < 3; ++d) {
                if (xi[d] < box_min[d]) box_min[d] = xi[d];
                if (xi[d] > box_max[d]) box_max[d] = xi[d];
            }
        }
        // Slight inflation so all particles fall strictly inside the bbox.
        double extent[3];
        double scale[3];
        constexpr std::uint64_t MORTON_RES = (1ull << 21) - 1;
        for (int d = 0; d < 3; ++d) {
            extent[d] = std::max(box_max[d] - box_min[d], 1e-12);
            scale[d] = static_cast<double>(MORTON_RES) / extent[d];
        }
        // 2. Build (morton_code, particle_index) pairs in serial.
        static std::vector<std::pair<std::uint64_t, std::uint32_t>> morton_pairs;
        if (morton_pairs.size() < N) morton_pairs.resize(N);
        for (std::uint32_t i = 0; i < N; ++i) {
            const double* xi = x + i * 3u;
            const std::uint32_t ix = static_cast<std::uint32_t>(
                std::min((xi[0] - box_min[0]) * scale[0],
                         static_cast<double>(MORTON_RES)));
            const std::uint32_t iy = static_cast<std::uint32_t>(
                std::min((xi[1] - box_min[1]) * scale[1],
                         static_cast<double>(MORTON_RES)));
            const std::uint32_t iz = static_cast<std::uint32_t>(
                std::min((xi[2] - box_min[2]) * scale[2],
                         static_cast<double>(MORTON_RES)));
            morton_pairs[i] = {rcp_morton_3d_encode(ix, iy, iz), i};
        }
#if RCP_HAVE_GNU_PARALLEL_SORT
        __gnu_parallel::sort(morton_pairs.begin(), morton_pairs.begin() + N);
#else
        std::sort(morton_pairs.begin(), morton_pairs.begin() + N);
#endif
        // 3. Extract sorted perm and sorted_morton arrays in serial.
        if (sorted_morton.size() < N) sorted_morton.resize(N);
        for (std::uint32_t i = 0; i < N; ++i) {
            perm[i] = morton_pairs[i].second;
            sorted_morton[i] = morton_pairs[i].first;
        }
        // 4. Recursive build. Top 3 bits at shift 60 define root's 8 octants.
        // Max Morton bit position is 62 (63 bits total: 21 per axis × 3).
        // For octree level 0 use bits [60, 62], level 1 uses [57, 59], etc.
        build_recursive_lbvh(x, D_arr, 0,
                             0, N, 0,
                             /*morton_shift=*/60u);
    }

    // Cycle 11 H38: octree build via Morton-bit-triplet octant boundaries.
    // perm[] is sorted by Morton code, so particles in the same octant at a
    // given depth form a contiguous range. Scan the range to find the 8
    // boundaries.
    void build_recursive_lbvh(const double* x, const double* D_arr,
                              std::uint32_t node_idx,
                              std::uint32_t begin, std::uint32_t end,
                              std::uint32_t depth,
                              std::uint32_t morton_shift) {
        OctNode& node = nodes[node_idx];
        depth_of_node[node_idx] = static_cast<std::uint8_t>(
            depth > 255u ? 255u : depth);
        const std::uint32_t count = end - begin;
        // Leaf: stop when count is small enough OR we've exhausted Morton
        // bits (extremely deep — shouldn't happen at realistic N).
        if (count <= KD_LEAF_SIZE || morton_shift < 3) {
            double lo[3] = {std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()};
            double hi[3] = {-std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()};
            double mx_D = 0.0;
            for (std::uint32_t k = begin; k < end; ++k) {
                const std::uint32_t p = perm[k];
                const double* xp = x + p * 3u;
                for (std::size_t d = 0; d < 3; ++d) {
                    const double v = xp[d];
                    if (v < lo[d]) lo[d] = v;
                    if (v > hi[d]) hi[d] = v;
                }
                const double Dp = D_arr[p];
                if (Dp > mx_D) mx_D = Dp;
                leaf_of_particle[p] = node_idx;
            }
            for (std::size_t d = 0; d < 3; ++d) {
                node.bbox_lo[d] = rcp_f_round_down(lo[d]);
                node.bbox_hi[d] = rcp_f_round_up(hi[d]);
            }
            node.max_D = rcp_f_round_up(mx_D);
            node.first = begin;
            node.count = count;
            for (int k = 0; k < 8; ++k) node.children[k] = KD_NULL;
            return;
        }
        // Internal: find octant boundaries by scanning Morton bit triplet.
        std::uint32_t ob[9];
        std::uint32_t pos = begin;
        for (int oct = 0; oct < 8; ++oct) {
            ob[oct] = pos;
            while (pos < end &&
                   ((sorted_morton[pos] >> morton_shift) & 0x7ull) ==
                       static_cast<std::uint64_t>(oct)) {
                ++pos;
            }
        }
        ob[8] = pos;
        // pos should equal end here (Morton sort guarantees all 8 octant IDs
        // 0..7 are accounted for in ascending order). If a few particles
        // remain (shouldn't happen with correctly-sorted Morton), they get
        // dropped — defensive: extend octant 7 to cover anything trailing.
        if (pos < end) ob[8] = end;
        // Allocate 8 child slots, recurse into non-empty octants.
        const std::uint32_t children_start =
            next_node.fetch_add(8, std::memory_order_relaxed);
        for (int k = 0; k < 8; ++k) {
            const std::uint32_t cb = ob[k];
            const std::uint32_t ce = ob[k + 1];
            if (cb == ce) {
                node.children[k] = KD_NULL;
                continue;
            }
            const std::uint32_t ci = children_start + k;
            node.children[k] = ci;
            parent_of_node[ci] = node_idx;
            build_recursive_lbvh(x, D_arr, ci, cb, ce, depth + 1,
                                 morton_shift - 3);
        }
        // Internal bbox/max_D = union of children.
        float lo[3] = {std::numeric_limits<float>::infinity(),
                       std::numeric_limits<float>::infinity(),
                       std::numeric_limits<float>::infinity()};
        float hi[3] = {-std::numeric_limits<float>::infinity(),
                       -std::numeric_limits<float>::infinity(),
                       -std::numeric_limits<float>::infinity()};
        float mx_D = 0.0f;
        for (int k = 0; k < 8; ++k) {
            const std::uint32_t ci = node.children[k];
            if (ci == KD_NULL) continue;
            const OctNode& c = nodes[ci];
            for (std::size_t d = 0; d < 3; ++d) {
                if (c.bbox_lo[d] < lo[d]) lo[d] = c.bbox_lo[d];
                if (c.bbox_hi[d] > hi[d]) hi[d] = c.bbox_hi[d];
            }
            if (c.max_D > mx_D) mx_D = c.max_D;
        }
        for (std::size_t d = 0; d < 3; ++d) {
            node.bbox_lo[d] = lo[d];
            node.bbox_hi[d] = hi[d];
        }
        node.max_D = mx_D;
        node.first = begin;
        node.count = 0;  // mark internal
    }

    static constexpr std::uint32_t OCT_PARALLEL_SPLIT_MIN = 32768;

    // Split perm[begin..end) into 8 octants via 3 nested median splits
    // (X, then Y within each X-half, then Z within each XY-quadrant).
    // Writes the boundary indices into oct_bound[0..8]; oct_bound[0]=begin,
    // oct_bound[8]=end, intermediate values are the 7 split points.
    // Cycle 10 H43-E2: when use_parallel is true (root level), use
    // __gnu_parallel::nth_element for the heavy first splits. This is the
    // safe alternative to OMP tasks: parallelizes the dominant serial cost
    // without the task-overhead pitfall that hung Phase E at N=500K.
    static void partition_into_octants(
        std::vector<std::uint32_t>& perm, const double* x,
        std::uint32_t begin, std::uint32_t end,
        std::uint32_t oct_bound[9],
        bool use_parallel) {
        const std::uint32_t mid_x = begin + (end - begin) / 2;
        auto by_x = [x](std::uint32_t a, std::uint32_t b) {
            return x[a * 3] < x[b * 3];
        };
        auto by_y = [x](std::uint32_t a, std::uint32_t b) {
            return x[a * 3 + 1] < x[b * 3 + 1];
        };
#if RCP_HAVE_GNU_PARALLEL_SORT
        if (use_parallel) {
            __gnu_parallel::nth_element(
                perm.begin() + begin, perm.begin() + mid_x,
                perm.begin() + end, by_x);
        } else {
            std::nth_element(
                perm.begin() + begin, perm.begin() + mid_x,
                perm.begin() + end, by_x);
        }
#else
        std::nth_element(
            perm.begin() + begin, perm.begin() + mid_x,
            perm.begin() + end, by_x);
#endif
        const std::uint32_t mid_y_lo = begin + (mid_x - begin) / 2;
#if RCP_HAVE_GNU_PARALLEL_SORT
        if (use_parallel) {
            __gnu_parallel::nth_element(
                perm.begin() + begin, perm.begin() + mid_y_lo,
                perm.begin() + mid_x, by_y);
        } else {
            std::nth_element(
                perm.begin() + begin, perm.begin() + mid_y_lo,
                perm.begin() + mid_x, by_y);
        }
#else
        std::nth_element(
            perm.begin() + begin, perm.begin() + mid_y_lo,
            perm.begin() + mid_x, by_y);
#endif
        const std::uint32_t mid_y_hi = mid_x + (end - mid_x) / 2;
#if RCP_HAVE_GNU_PARALLEL_SORT
        if (use_parallel) {
            __gnu_parallel::nth_element(
                perm.begin() + mid_x, perm.begin() + mid_y_hi,
                perm.begin() + end, by_y);
        } else {
            std::nth_element(
                perm.begin() + mid_x, perm.begin() + mid_y_hi,
                perm.begin() + end, by_y);
        }
#else
        std::nth_element(
            perm.begin() + mid_x, perm.begin() + mid_y_hi,
            perm.begin() + end, by_y);
#endif
        // 4 quadrants now: [begin, mid_y_lo), [mid_y_lo, mid_x),
        //                   [mid_x, mid_y_hi),  [mid_y_hi, end).
        // Split each by Z.
        const std::uint32_t mid_z_q0 = begin + (mid_y_lo - begin) / 2;
        const std::uint32_t mid_z_q1 = mid_y_lo + (mid_x - mid_y_lo) / 2;
        const std::uint32_t mid_z_q2 = mid_x + (mid_y_hi - mid_x) / 2;
        const std::uint32_t mid_z_q3 = mid_y_hi + (end - mid_y_hi) / 2;
        auto by_z = [x](std::uint32_t a, std::uint32_t b) {
            return x[a * 3 + 2] < x[b * 3 + 2];
        };
        std::nth_element(perm.begin() + begin, perm.begin() + mid_z_q0,
                         perm.begin() + mid_y_lo, by_z);
        std::nth_element(perm.begin() + mid_y_lo, perm.begin() + mid_z_q1,
                         perm.begin() + mid_x, by_z);
        std::nth_element(perm.begin() + mid_x, perm.begin() + mid_z_q2,
                         perm.begin() + mid_y_hi, by_z);
        std::nth_element(perm.begin() + mid_y_hi, perm.begin() + mid_z_q3,
                         perm.begin() + end, by_z);
        oct_bound[0] = begin;
        oct_bound[1] = mid_z_q0;
        oct_bound[2] = mid_y_lo;
        oct_bound[3] = mid_z_q1;
        oct_bound[4] = mid_x;
        oct_bound[5] = mid_z_q2;
        oct_bound[6] = mid_y_hi;
        oct_bound[7] = mid_z_q3;
        oct_bound[8] = end;
    }

    void build_recursive_oct(const double* x, const double* D_arr,
                             std::uint32_t node_idx,
                             std::uint32_t begin, std::uint32_t end,
                             std::uint32_t depth) {
        OctNode& node = nodes[node_idx];
        depth_of_node[node_idx] = static_cast<std::uint8_t>(
            depth > 255u ? 255u : depth);
        const std::uint32_t count = end - begin;
        if (count <= KD_LEAF_SIZE) {
            // Leaf: compute tight bbox + max_D from particles.
            double lo[3] = {std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()};
            double hi[3] = {-std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()};
            double mx_D = 0.0;
            for (std::uint32_t k = begin; k < end; ++k) {
                const std::uint32_t p = perm[k];
                const double* xp = x + p * 3u;
                for (std::size_t d = 0; d < 3; ++d) {
                    const double v = xp[d];
                    if (v < lo[d]) lo[d] = v;
                    if (v > hi[d]) hi[d] = v;
                }
                const double Dp = D_arr[p];
                if (Dp > mx_D) mx_D = Dp;
                leaf_of_particle[p] = node_idx;
            }
            for (std::size_t d = 0; d < 3; ++d) {
                node.bbox_lo[d] = rcp_f_round_down(lo[d]);
                node.bbox_hi[d] = rcp_f_round_up(hi[d]);
            }
            node.max_D = rcp_f_round_up(mx_D);
            node.first = begin;
            node.count = count;
            for (int k = 0; k < 8; ++k) node.children[k] = KD_NULL;
            return;
        }
        // Internal: partition into 8 octants and recurse.
        std::uint32_t ob[9];
        // H43-E2 root-parallel tested + reverted (caused hangs at N=500K
        // identical to H43-E). Octree's static state may be incompatible
        // with __gnu_parallel's internal OpenMP region. Fall back to fully
        // serial build; future parallelism via LBVH (Morton-sort) instead.
        partition_into_octants(perm, x, begin, end, ob,
                               /*use_parallel=*/false);
        const std::uint32_t children_start =
            next_node.fetch_add(8, std::memory_order_relaxed);
        // Recurse into each non-empty octant. Empty octants get marked as
        // KD_NULL so range_query can skip them quickly.
        // Cycle 10 H43-E: spawn OMP tasks for child subtrees ≥ threshold so
        // parallel build matches binary tree's H19+ pattern.
        if (count >= OCT_PARALLEL_SPLIT_MIN) {
            for (int k = 0; k < 8; ++k) {
                const std::uint32_t child_begin = ob[k];
                const std::uint32_t child_end = ob[k + 1];
                if (child_begin == child_end) {
                    node.children[k] = KD_NULL;
                    continue;
                }
                const std::uint32_t child_idx = children_start + k;
                node.children[k] = child_idx;
                parent_of_node[child_idx] = node_idx;
                #pragma omp task firstprivate(child_idx, child_begin, child_end) shared(x, D_arr)
                build_recursive_oct(x, D_arr, child_idx,
                                    child_begin, child_end, depth + 1);
            }
        } else {
            for (int k = 0; k < 8; ++k) {
                const std::uint32_t child_begin = ob[k];
                const std::uint32_t child_end = ob[k + 1];
                if (child_begin == child_end) {
                    node.children[k] = KD_NULL;
                    continue;
                }
                const std::uint32_t child_idx = children_start + k;
                node.children[k] = child_idx;
                parent_of_node[child_idx] = node_idx;
                build_recursive_oct(x, D_arr, child_idx,
                                    child_begin, child_end, depth + 1);
            }
        }
        // If we spawned tasks, wait for them before computing internal bbox.
        if (count >= OCT_PARALLEL_SPLIT_MIN) {
            #pragma omp taskwait
        }
        // Compute internal bbox as union of present children.
        float lo[3] = {std::numeric_limits<float>::infinity(),
                       std::numeric_limits<float>::infinity(),
                       std::numeric_limits<float>::infinity()};
        float hi[3] = {-std::numeric_limits<float>::infinity(),
                       -std::numeric_limits<float>::infinity(),
                       -std::numeric_limits<float>::infinity()};
        float mx_D = 0.0f;
        for (int k = 0; k < 8; ++k) {
            const std::uint32_t ci = node.children[k];
            if (ci == KD_NULL) continue;
            const OctNode& c = nodes[ci];
            for (std::size_t d = 0; d < 3; ++d) {
                if (c.bbox_lo[d] < lo[d]) lo[d] = c.bbox_lo[d];
                if (c.bbox_hi[d] > hi[d]) hi[d] = c.bbox_hi[d];
            }
            if (c.max_D > mx_D) mx_D = c.max_D;
        }
        for (std::size_t d = 0; d < 3; ++d) {
            node.bbox_lo[d] = lo[d];
            node.bbox_hi[d] = hi[d];
        }
        node.max_D = mx_D;
        node.first = begin;
        node.count = 0;  // mark internal
    }

    // Recompute a leaf's bbox + max_D after particles in it moved.
    // Mirror of KdTree::recompute_leaf_bbox.
    void recompute_leaf_bbox(std::uint32_t leaf_idx,
                             const double* x, const double* D_arr) {
        OctNode& L = nodes[leaf_idx];
        double lo[3] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[3] = {-std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
        double mx_D = 0.0;
        const std::uint32_t end = L.first + L.count;
        for (std::uint32_t k = L.first; k < end; ++k) {
            const std::uint32_t p = perm[k];
            const double* xp = x + p * 3u;
            for (std::size_t d = 0; d < 3; ++d) {
                if (xp[d] < lo[d]) lo[d] = xp[d];
                if (xp[d] > hi[d]) hi[d] = xp[d];
            }
            const double Dk = D_arr[p];
            if (Dk > mx_D) mx_D = Dk;
        }
        for (std::size_t d = 0; d < 3; ++d) {
            L.bbox_lo[d] = rcp_f_round_down(lo[d]);
            L.bbox_hi[d] = rcp_f_round_up(hi[d]);
        }
        L.max_D = rcp_f_round_up(mx_D);
    }

    // 8-way scalar range query. SIMD comes in Phase D.
    void range_query(const double* xq, double Di, double t_shell,
                     const double* x, const double* D_arr,
                     std::vector<std::uint32_t>& out) const {
        if (nodes.empty()) return;
        constexpr std::size_t MAX_STACK = 256;
        std::uint32_t stack[MAX_STACK];
        int top = 0;
        stack[top++] = 0;
        const double xq0 = xq[0], xq1 = xq[1], xq2 = xq[2];
        while (top > 0) {
            const std::uint32_t node_idx = stack[--top];
            const OctNode& n = nodes[node_idx];
            // Bbox prune for this node.
            const double diff0 = std::max(0.0,
                std::max((double)n.bbox_lo[0] - xq0,
                         xq0 - (double)n.bbox_hi[0]));
            const double diff1 = std::max(0.0,
                std::max((double)n.bbox_lo[1] - xq1,
                         xq1 - (double)n.bbox_hi[1]));
            const double diff2 = std::max(0.0,
                std::max((double)n.bbox_lo[2] - xq2,
                         xq2 - (double)n.bbox_hi[2]));
            const double dbb2 = diff0*diff0 + diff1*diff1 + diff2*diff2;
            const double r_node = (Di + n.max_D) * 0.5 + t_shell;
            if (dbb2 > r_node * r_node) continue;
            if (n.count > 0) {
                // Leaf — scan particles.
                const std::uint32_t cnt = n.count;
                const std::uint32_t first = n.first;
                for (std::uint32_t k = 0; k < cnt; ++k) {
                    const std::uint32_t p = perm[first + k];
                    const double* xp = x + p * 3u;
                    const double dd0 = xp[0] - xq0;
                    const double dd1 = xp[1] - xq1;
                    const double dd2 = xp[2] - xq2;
                    const double d2 = dd0*dd0 + dd1*dd1 + dd2*dd2;
                    const double r_p = (Di + D_arr[p]) * 0.5 + t_shell;
                    if (d2 <= r_p * r_p) out.push_back(p);
                }
                continue;
            }
            // Internal — push non-null children (scalar; Phase D SIMDs this).
            for (int k = 0; k < 8; ++k) {
                if (n.children[k] != KD_NULL) {
                    stack[top++] = n.children[k];
                }
            }
        }
    }
};

// Per-axis PBC image offsets for a query at xi[d] with search radius r against
// a box of size box[d]. Returns the offsets {-1, 0, +1} that need to be
// considered (i.e. shifted query sphere overlaps [0, box[d]) ).
struct ImageSet {
    std::int8_t offsets[KD_MAX_NDIM][3];
    std::uint8_t counts[KD_MAX_NDIM];
    std::uint32_t total_images = 1;
};

static ImageSet compute_image_set(const double* xi, double r,
                                  const double* box, std::size_t Ndim) {
    // Query at xq = xi + m·box, looking for particles whose unwrapped
    // position is within r of xq. The relevant offset m for axis d is:
    //   - m = 0 always (look for nearby in same cell)
    //   - m = +1 when xi[d] < r — to find particles that wrap around from
    //     the FAR (high) side of the box, because their MIC distance is
    //     achieved by treating them at x[j] - box (i.e., "before" the box).
    //     Equivalently, query at xi + box hits those particles' actual
    //     positions in [xi + box - r, xi + box + r] ∩ [0, box[d]).
    //   - m = -1 when xi[d] > box - r — symmetric: find particles that
    //     wrap around from the NEAR (low) side.
    ImageSet s;
    s.total_images = 1;
    for (std::size_t d = 0; d < Ndim; ++d) {
        std::uint8_t c = 0;
        s.offsets[d][c++] = 0;
        if (xi[d] < r) s.offsets[d][c++] = +1;
        if (xi[d] > box[d] - r) s.offsets[d][c++] = -1;
        s.counts[d] = c;
        s.total_images *= c;
    }
    return s;
}

}  // namespace

void get_pairs_kdtree(
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
    std::uint32_t* out_max_observed,
    std::uint32_t* out_min_observed)
{
    if (N == 0 || Ndim == 0) {
        if (out_max_observed) *out_max_observed = 0;
        if (out_min_observed) *out_min_observed = 0;
        (void)walls;
        return;
    }
    // kdtree Cycle 4: probe bracket — master-thread wall covers all phases.
    const bool kd_probe_on = rcp_cycle_probe_enabled();
    const std::uint64_t kd_tsc_fn_start = kd_probe_on ? __rdtsc() : 0;
    const std::uint64_t kd_tsc_setup_start = kd_tsc_fn_start;

    // ----- D-stats with index caching (Cycle 6 H4) --------------------------
    // D[i] = D0[i] * kappa with scalar kappa, so the indices of max/min/
    // median are invariant across calls under uniform rescale. Cache the
    // indices and do O(1) lookups instead of 3×O(N) reductions.
    //
    // Spatial reorder (driver Cycle 10) permutes the D[] array → invalidates
    // cached indices. Detect by tracking two sentinel D values: under
    // uniform scaling, D[0]/cached_D0 == D[N/2]/cached_Dmid. After a reorder
    // these ratios diverge dramatically. ULP-tolerant equality (1e-9 rel)
    // accommodates FP rounding while catching any real permutation.
    double Dmax, Dmin, median;
    {
        static std::size_t cached_max_idx = 0;
        static std::size_t cached_min_idx = 0;
        static std::size_t cached_med_idx = 0;
        static std::size_t cached_N = 0;
        static double cached_D0 = 0.0, cached_Dmid = 0.0;
        static bool cached_valid = false;
        bool hit = false;
        if (cached_valid && cached_N == N && cached_D0 != 0.0 &&
            cached_Dmid != 0.0) {
            const double s0   = D[0]   / cached_D0;
            const double smid = D[N/2] / cached_Dmid;
            const double tol  = 1e-9 * std::max(std::abs(s0), 1.0);
            if (std::abs(s0 - smid) <= tol) {
                Dmax   = D[cached_max_idx];
                Dmin   = D[cached_min_idx];
                median = D[cached_med_idx];
                cached_D0   = D[0];
                cached_Dmid = D[N/2];
                hit = true;
            }
        }
        if (!hit) {
            auto max_it = std::max_element(D.begin(), D.end());
            auto min_it = std::min_element(D.begin(), D.end());
            cached_max_idx = static_cast<std::size_t>(max_it - D.begin());
            cached_min_idx = static_cast<std::size_t>(min_it - D.begin());
            Dmax = *max_it;
            Dmin = *min_it;
            std::vector<std::size_t> idx_sorted(N);
            std::iota(idx_sorted.begin(), idx_sorted.end(),
                      static_cast<std::size_t>(0));
            std::nth_element(idx_sorted.begin(),
                             idx_sorted.begin() + N / 2, idx_sorted.end(),
                             [&D](std::size_t a, std::size_t b) {
                                 return D[a] < D[b];
                             });
            cached_med_idx = idx_sorted[N / 2];
            median = D[cached_med_idx];
            cached_N    = N;
            cached_D0   = D[0];
            cached_Dmid = D[N/2];
            cached_valid = true;
        }
    }
    double t_shell = std::min(std::max(median * 1.02, Dmin * 2.25), Dmax) * 0.5;
    // Cycle 14 H53: tunable scale knob to sweep shell thickness. Default 1.0.
    // RCP_T_SHELL_SCALE=0.5 halves the shell margin → tighter lists → less
    // force.loop & tree work, but risks missed interactions if drift exceeds
    // the new margin between K-batched get_pairs calls.
    {
        static const double t_shell_scale = []() {
            const char* s = std::getenv("RCP_T_SHELL_SCALE");
            if (!s || !s[0]) return 1.0;
            double v = std::atof(s);
            return (v > 0.0 && v < 10.0) ? v : 1.0;
        }();
        t_shell *= t_shell_scale;
    }

    // Precompute inv_box for the MIC math — same as get_pairs_nd_3 so that
    // borderline pairs (d ≈ r_c) classify identically.
    double inv_box_kd[KD_MAX_NDIM] = {0.0};
    for (std::size_t d = 0; d < Ndim; ++d) inv_box_kd[d] = 1.0 / box[d];

    // ----- Per-particle locks for cross-writes (mirror sort+walk) -----------
    static std::vector<omp_lock_t> j_locks_kd;
    static std::size_t j_locks_kd_N = 0;
    if (j_locks_kd_N != N) {
        for (std::size_t k = 0; k < j_locks_kd_N; ++k)
            omp_destroy_lock(&j_locks_kd[k]);
        j_locks_kd.resize(N);
        for (std::size_t k = 0; k < N; ++k) omp_init_lock(&j_locks_kd[k]);
        j_locks_kd_N = N;
    }

    // kdtree Cycle 6 (H1): dirty-bit per particle. A list needs sorting only
    // when it was touched this call — either own-write (refresh[i]=1 rebuild)
    // or cross-write (some refreshed j appended to slot_i). Initialize from
    // refresh[]; cross-writes set dirty_kd[j]=1 alongside the slot append.
    //
    // Cycle 10 H29: also build the sparse refresh-list in the same pass so
    // downstream sites (patch, probe, main parallel-for) iterate R ≪ N
    // entries instead of scanning all N. In mid regime R/N ≈ 0.5% — the
    // O(N) scans were dominating per-call fixed cost.
    // Cycle 10 H31: also maintain kd_dirty_list (refresh-dirty + cross-write
    // dirty). The sort pass then iterates the small dirty list instead of
    // scanning all N. Initially populated with refreshed particles; cross-
    // write recipients append themselves during the serial merge below.
    // Cycle 12 H47 (parallel dirty_kd init + refresh-list build) tested at
    // N=100K and N=500K, both showed net +5-15s regression vs serial. The
    // OMP region setup cost per call (~10-30us) at our 10K call rate
    // exceeded the parallel savings on this small O(N) loop. REVERTED to
    // serial — simpler and faster at our workload's call rate.
    static std::vector<std::uint8_t> dirty_kd;
    static std::vector<std::uint32_t> kd_refresh_list;
    static std::vector<std::uint32_t> kd_dirty_list;
    if (dirty_kd.size() != N) dirty_kd.resize(N);
    kd_refresh_list.clear();
    kd_dirty_list.clear();
    if (kd_refresh_list.capacity() < N) kd_refresh_list.reserve(N);
    if (kd_dirty_list.capacity() < N) kd_dirty_list.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        if (refresh[k] == 1) {
            dirty_kd[k] = 1;
            kd_refresh_list.push_back(static_cast<std::uint32_t>(k));
            kd_dirty_list.push_back(static_cast<std::uint32_t>(k));
        } else {
            dirty_kd[k] = 0;
        }
    }
    const std::size_t kd_R = kd_refresh_list.size();

    if (kd_probe_on) {
        g_kd_setup_cycles.fetch_add(__rdtsc() - kd_tsc_setup_start,
                                    std::memory_order_relaxed);
    }

    // ----- Build kd-tree or patch (Cycle 7 H10: adaptive) -------------------
    // Strategy from instrumented 25K-step probe data: early regime (high
    // refresh count, ~all leaves touched) costs more in queries than build
    // and has zero patching benefit; mid/late regime (low refresh, ~10% of
    // leaves touched, queries cheap) inverts to build-dominated and gains
    // big from patching. Use refresh count as the regime signal: above
    // threshold → full rebuild and tight cutoffs; below → patch refreshed
    // leaves and inflate query radius by Dmin/2 to cover stale leaves.
    // K_MAX_PATCHES forces a periodic full rebuild for safety against
    // accumulating staleness in never-patched leaves.
    // Cycle 12 H45: a-priori model says original 1/20 (5%) is well-tuned.
    // Patch path adds Dmin/2 query-radius inflation, increasing per-query
    // cost by ~30%. For refresh > 5%, this inflation × refresh-count exceeds
    // the build-cost savings even with LBVH-fast rebuild. Leave at 5%.
    constexpr std::size_t KD_PATCH_REFRESH_THRESHOLD_FRAC = 20; // 1/20 = 5%
    constexpr std::size_t KD_MAX_PATCHES_BETWEEN_REBUILDS = 256;
    static KdTree tree;
    static OctTree oct_tree;                                  // Cycle 10 H43
    static std::size_t g_kd_patches_since_rebuild = 0;
    static std::size_t g_kd_tree_N = 0;
    const bool use_octree = rcp_use_octree();                 // Cycle 10 H43
    // Cycle 10 H29: refresh count taken from the sparse refresh-list (built
    // above), eliminating an extra O(N) scan.
    const std::size_t refresh_count_total = kd_R;
    // Cycle 10 H43: shared rebuild trigger for both tree variants. Phase F:
    // patch path for octree now supported.
    const bool need_full_rebuild =
        (g_kd_tree_N != N) ||
        (g_kd_patches_since_rebuild >= KD_MAX_PATCHES_BETWEEN_REBUILDS) ||
        (refresh_count_total * KD_PATCH_REFRESH_THRESHOLD_FRAC > N);

    const std::uint64_t kd_tsc_build_start = kd_probe_on ? __rdtsc() : 0;
    if (need_full_rebuild) {
        if (use_octree) {
            oct_tree.build(x, D.data(), N, Ndim);
        } else {
            tree.build(x, D.data(), N, Ndim);
        }
        g_kd_tree_N = N;
        g_kd_patches_since_rebuild = 0;
    } else if (use_octree) {
        // Cycle 10 H43 Phase F: octree patch path adapted from H30. Same
        // structure: collect distinct touched leaves + recompute bboxes,
        // then bottom-up by depth recompute ancestor bboxes (8-way union).
        const std::size_t node_count = oct_tree.nodes.size();
        static std::vector<std::uint8_t> g_oct_node_touched;
        if (g_oct_node_touched.size() < node_count)
            g_oct_node_touched.resize(node_count);
        std::fill(g_oct_node_touched.begin(),
                  g_oct_node_touched.begin() + node_count,
                  std::uint8_t{0});
        constexpr int OCT_MAX_DEPTH = 16;
        static std::vector<std::uint32_t> g_oct_depth_bucket[OCT_MAX_DEPTH];
        int max_depth_seen = -1;
        for (std::size_t r = 0; r < kd_R; ++r) {
            const std::size_t i = kd_refresh_list[r];
            const std::uint32_t leaf = oct_tree.leaf_of_particle[i];
            if (g_oct_node_touched[leaf]) continue;
            g_oct_node_touched[leaf] = 1;
            oct_tree.recompute_leaf_bbox(leaf, x, D.data());
            std::uint32_t curr = oct_tree.parent_of_node[leaf];
            while (curr != KD_NULL && !g_oct_node_touched[curr]) {
                g_oct_node_touched[curr] = 1;
                const int d = static_cast<int>(oct_tree.depth_of_node[curr]);
                if (d >= 0 && d < OCT_MAX_DEPTH) {
                    g_oct_depth_bucket[d].push_back(curr);
                    if (d > max_depth_seen) max_depth_seen = d;
                }
                curr = oct_tree.parent_of_node[curr];
            }
        }
        for (int d = max_depth_seen; d >= 0; --d) {
            auto& bucket = g_oct_depth_bucket[d];
            for (std::uint32_t idx : bucket) {
                OctNode& n = oct_tree.nodes[idx];
                float lo[3] = {std::numeric_limits<float>::infinity(),
                               std::numeric_limits<float>::infinity(),
                               std::numeric_limits<float>::infinity()};
                float hi[3] = {-std::numeric_limits<float>::infinity(),
                               -std::numeric_limits<float>::infinity(),
                               -std::numeric_limits<float>::infinity()};
                float mx_D = 0.0f;
                for (int k = 0; k < 8; ++k) {
                    const std::uint32_t ci = n.children[k];
                    if (ci == KD_NULL) continue;
                    const OctNode& c = oct_tree.nodes[ci];
                    for (int dd = 0; dd < 3; ++dd) {
                        if (c.bbox_lo[dd] < lo[dd]) lo[dd] = c.bbox_lo[dd];
                        if (c.bbox_hi[dd] > hi[dd]) hi[dd] = c.bbox_hi[dd];
                    }
                    if (c.max_D > mx_D) mx_D = c.max_D;
                }
                for (int dd = 0; dd < 3; ++dd) {
                    n.bbox_lo[dd] = lo[dd];
                    n.bbox_hi[dd] = hi[dd];
                }
                n.max_D = mx_D;
            }
            bucket.clear();
        }
        ++g_kd_patches_since_rebuild;
    } else {
        // Patch each refreshed particle's leaf, then propagate ancestor
        // bboxes up. Cycle 10 H30: previous version called
        // tree.propagate_up_from once per touched leaf, redundantly
        // recomputing each ancestor up to (# touched leaves below it) times.
        // Replace with a two-phase: (A) collect distinct touched leaves and
        // recompute their bboxes; (B) BFS up by descending depth, each
        // ancestor visited exactly once. Cost goes from
        // O(L × depth) → O(distinct ancestors) ≤ O(L × log2(T/L)).
        const std::size_t node_count = tree.nodes.size();
        static std::vector<std::uint8_t> g_kd_node_touched;
        if (g_kd_node_touched.size() < node_count)
            g_kd_node_touched.resize(node_count);
        std::fill(g_kd_node_touched.begin(),
                  g_kd_node_touched.begin() + node_count,
                  std::uint8_t{0});
        // Per-depth buckets of node indices to recompute. uint8 depth caps at
        // 255 levels which is far beyond any realistic tree.
        constexpr int H30_MAX_DEPTH = 64;
        static std::vector<std::uint32_t> g_kd_depth_bucket[H30_MAX_DEPTH];
        int max_depth_seen = -1;
        // Phase A: distinct touched leaves + recompute their bboxes,
        // record their ancestors in depth buckets.
        for (std::size_t r = 0; r < kd_R; ++r) {
            const std::size_t i = kd_refresh_list[r];
            const std::uint32_t leaf = tree.leaf_of_particle[i];
            if (g_kd_node_touched[leaf]) continue;
            g_kd_node_touched[leaf] = 1;
            tree.recompute_leaf_bbox(leaf, x, D.data(), Ndim);
            // Walk up adding each unseen ancestor to its depth bucket.
            std::uint32_t curr = tree.parent_of_node[leaf];
            while (curr != KD_NULL && !g_kd_node_touched[curr]) {
                g_kd_node_touched[curr] = 1;
                const int d = static_cast<int>(tree.depth_of_node[curr]);
                if (d >= 0 && d < H30_MAX_DEPTH) {
                    g_kd_depth_bucket[d].push_back(curr);
                    if (d > max_depth_seen) max_depth_seen = d;
                }
                curr = tree.parent_of_node[curr];
            }
        }
        // Phase B: process by DECREASING depth so each ancestor sees its
        // children's already-current bboxes (whether touched or untouched).
        for (int d = max_depth_seen; d >= 0; --d) {
            auto& bucket = g_kd_depth_bucket[d];
            for (std::uint32_t idx : bucket) {
                KdNode& n = tree.nodes[idx];
                const KdNode& lc = tree.nodes[n.children[0]];
                const KdNode& rc = tree.nodes[n.children[1]];
                for (std::size_t dd = 0; dd < Ndim; ++dd) {
                    n.bbox_lo[dd] = std::min(lc.bbox_lo[dd], rc.bbox_lo[dd]);
                    n.bbox_hi[dd] = std::max(lc.bbox_hi[dd], rc.bbox_hi[dd]);
                }
                n.max_D = std::max(lc.max_D, rc.max_D);
            }
            bucket.clear();
        }
        ++g_kd_patches_since_rebuild;
    }
    const std::uint64_t kd_build_cycles_this_call =
        kd_probe_on ? (__rdtsc() - kd_tsc_build_start) : 0;
    if (kd_probe_on) {
        g_kd_build_cycles.fetch_add(kd_build_cycles_this_call,
                                    std::memory_order_relaxed);
    }
    // Cycle 7 H10: inflate query radius by Dmin/2 only on patch calls. On
    // rebuild call, tree's bboxes are tight — no inflation needed (preserves
    // H3's tight per-pair cutoff for expensive early-regime queries).
    // Cycle 13 b6.A Phase 2 (scaling inflation with K_BATCH) tested + reverted:
    // bigger inflation grew query radius → more candidate evaluation per
    // query → wall went UP. And at K=4 phi-final dropped from 0.846 → 0.823,
    // suggesting the extra candidates perturbed convergence. Phase 1's
    // fixed Dmin/2 is empirically the best.
    const double kd_radius_inflation =
        need_full_rebuild ? 0.0 : (Dmin * 0.5);

    // Cycle 7 measurement: count refresh + leaves touched this call so we
    // can characterize the patching opportunity. Cheap O(N) scan; only when
    // probes on.
    std::uint64_t kd_refresh_count = 0;
    std::uint64_t kd_leaves_touched = 0;
    if (kd_probe_on) {
        // Reuse the static dirty_kd flag arrangement; count refresh and
        // distinct leaves via a small bitmap of leaf indices.
        // Cycle 10 H29: iterate sparse refresh-list instead of all N.
        // Cycle 10 H43: query the active tree variant for its nodes.size()
        // and per-particle leaf assignment.
        static std::vector<std::uint8_t> kd_leaf_seen_for_measure;
        const std::size_t node_count = use_octree
            ? oct_tree.nodes.size() : tree.nodes.size();
        if (kd_leaf_seen_for_measure.size() < node_count)
            kd_leaf_seen_for_measure.resize(node_count);
        std::fill(kd_leaf_seen_for_measure.begin(),
                  kd_leaf_seen_for_measure.begin() + node_count,
                  std::uint8_t{0});
        kd_refresh_count = static_cast<std::uint64_t>(kd_R);
        for (std::size_t r = 0; r < kd_R; ++r) {
            const std::size_t i = kd_refresh_list[r];
            const std::uint32_t leaf = use_octree
                ? oct_tree.leaf_of_particle[i]
                : tree.leaf_of_particle[i];
            if (!kd_leaf_seen_for_measure[leaf]) {
                kd_leaf_seen_for_measure[leaf] = 1;
                ++kd_leaves_touched;
            }
        }
        g_kd_refresh_hist[rcp_kd_log2_bin(kd_refresh_count)]
            .fetch_add(1, std::memory_order_relaxed);
        g_kd_leaves_touched_hist[rcp_kd_log2_bin(kd_leaves_touched)]
            .fetch_add(1, std::memory_order_relaxed);
        const int regime = rcp_kd_regime_for_step(g_shadow_current_step);
        g_kd_regime_calls[regime].fetch_add(1, std::memory_order_relaxed);
        g_kd_regime_refresh_sum[regime].fetch_add(kd_refresh_count,
                                                  std::memory_order_relaxed);
        g_kd_regime_leaves_touched_sum[regime].fetch_add(kd_leaves_touched,
                                                  std::memory_order_relaxed);
        g_kd_regime_build_cyc[regime].fetch_add(kd_build_cycles_this_call,
                                                  std::memory_order_relaxed);
    }

    bool warningIssued = false;

    // Cycle 7 H11: per-thread deferred cross-write buffer. Each thread
    // appends (j, i) pairs locally instead of locking slot_j. After the
    // parallel region, a single-threaded merge writes them to slot_j. The
    // existing sort+unique pass on dirty lists handles dedup. Eliminates
    // lock contention on big-particle slot_j in the early regime.
    const int kd_max_threads = omp_get_max_threads();
    static std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>
        g_kd_xwrite_buffers;
    if (static_cast<int>(g_kd_xwrite_buffers.size()) < kd_max_threads)
        g_kd_xwrite_buffers.resize(kd_max_threads);

    // ----- For each refreshed i, query tree, then enforce asymmetric write --
    const std::uint64_t kd_tsc_qregion_start = kd_probe_on ? __rdtsc() : 0;
    #pragma omp parallel shared(warningIssued)
    {
        std::vector<std::uint32_t> candidates;
        candidates.reserve(256);
        double xq[KD_MAX_NDIM];

        // Cycle 7 H11: per-thread cross-write buffer.
        auto& xwrite_buf = g_kd_xwrite_buffers[omp_get_thread_num()];
        xwrite_buf.clear();
        xwrite_buf.reserve(2048);

        // kdtree Cycle 4: thread-local probe accumulators.
        std::uint64_t tls_kd_image = 0;
        std::uint64_t tls_kd_query = 0;
        std::uint64_t tls_kd_mic   = 0;
        std::uint64_t tls_kd_write = 0;
        std::uint64_t tls_kd_parts = 0;
        std::uint64_t tls_kd_cand  = 0;
        std::uint64_t tls_kd_acc   = 0;
        std::uint64_t tls_kd_total = 0;
        const std::uint64_t tls_kd_tsc_thread_start = kd_probe_on ? __rdtsc() : 0;

        // Cycle 10 H29: drive the parallel-for over the sparse refresh-list
        // (kd_R entries) instead of all N. Every iteration now does real work
        // — no more `if (refresh[i]!=1) continue` skips eating scheduler
        // chunks. In mid regime (R ≪ N) this collapses N/32 = 15.6K chunk
        // grabs down to R/32 ≈ 80-300, eliminating most per-call fixed cost.
        #pragma omp for schedule(dynamic, 32)
        for (std::size_t r = 0; r < kd_R; ++r) {
            const std::size_t i = kd_refresh_list[r];
            if (kd_probe_on) ++tls_kd_parts;
            std::uint32_t* slot_i = pairs_data + pair_offsets[i];
            slot_i[0] = 0;
            const std::uint32_t cap_i = max_neighbors_per_particle[i];

            const double Di = D[i];
            const double r_c_max = (Di + 1.01 * Dmax) * 0.5 + t_shell;
            // Cycle 12 H51: r2_max was computed but never used. Removed.

            // Enumerate PBC images and run a range query per image.
            const std::uint64_t tsc_img0 = kd_probe_on ? __rdtsc() : 0;
            const double* xi = x + i * Ndim;
            const ImageSet img = compute_image_set(xi, r_c_max, box.data(), Ndim);
            if (kd_probe_on) tls_kd_image += __rdtsc() - tsc_img0;

            // Iterate the cartesian product of per-axis image offsets.
            // We unroll up to Ndim=8 here in a generic way using indices.
            std::size_t idx[KD_MAX_NDIM] = {0};
            for (std::uint32_t img_k = 0; img_k < img.total_images; ++img_k) {
                // Build xq for this image
                const std::uint64_t tsc_img1 = kd_probe_on ? __rdtsc() : 0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    const std::int8_t off = img.offsets[d][idx[d]];
                    xq[d] = xi[d] + static_cast<double>(off) * box[d];
                }
                if (kd_probe_on) tls_kd_image += __rdtsc() - tsc_img1;

                candidates.clear();
                const std::uint64_t tsc_q0 = kd_probe_on ? __rdtsc() : 0;
                // Cycle 6 H3: pass Di + (t_shell + inflation) so the tree
                // uses per-subtree max_D for tight pruning and per-pair
                // leaf cutoffs. The inflation handles drift of unrefreshed
                // particles in unpatched leaves; the outer MIC test below
                // re-tests against the strict cutoff using current MIC d.
                // Cycle 10 H43: dispatch to octree variant if active.
                if (use_octree) {
                    oct_tree.range_query(xq, Di,
                                         t_shell + kd_radius_inflation,
                                         x, D.data(), candidates);
                } else {
                    tree.range_query(xq, Di, t_shell + kd_radius_inflation,
                                     x, D.data(), Ndim, candidates);
                }
                if (kd_probe_on) {
                    tls_kd_query += __rdtsc() - tsc_q0;
                    tls_kd_cand  += candidates.size();
                }

                // For each candidate, apply inclusion + asymmetric storage.
                for (std::uint32_t jj : candidates) {
                    const std::size_t j = jj;
                    if (j == i) continue;  // self
                    // PBC-MIC distance check (kept after H3 because tree's
                    // leaf scan is on the post-image Euclidean distance for
                    // the queried image; the canonical MIC distance can be
                    // shorter when wrap is via a different image. Removing
                    // it costs more than it saves due to accept-write
                    // contention on the extra candidates.)
                    const std::uint64_t tsc_m0 = kd_probe_on ? __rdtsc() : 0;
                    double d2;
                    if (Ndim == 3) {
                        const double* xj = x + j * 3u;
                        double dd0 = xj[0] - xi[0];
                        double dd1 = xj[1] - xi[1];
                        double dd2 = xj[2] - xi[2];
                        dd0 -= std::nearbyint(dd0 * inv_box_kd[0]) * box[0];
                        dd1 -= std::nearbyint(dd1 * inv_box_kd[1]) * box[1];
                        dd2 -= std::nearbyint(dd2 * inv_box_kd[2]) * box[2];
                        d2 = dd0*dd0 + dd1*dd1 + dd2*dd2;
                    } else {
                        d2 = 0.0;
                        for (std::size_t d = 0; d < Ndim; ++d) {
                            double dd = x[j * Ndim + d] - xi[d];
                            dd -= std::nearbyint(dd * inv_box_kd[d]) * box[d];
                            d2 += dd * dd;
                        }
                    }
                    const double r_c = (Di + D[j]) * 0.5 + t_shell;
                    const bool reject = (d2 >= r_c * r_c);
                    if (kd_probe_on) tls_kd_mic += __rdtsc() - tsc_m0;
                    if (reject) continue;
                    if (kd_probe_on) ++tls_kd_acc;

                    const std::uint64_t tsc_w0 = kd_probe_on ? __rdtsc() : 0;
                    if (j > i) {
                        // Own-write to slot_i
                        std::uint32_t& cnt = slot_i[0];
                        // Dedup: check existing (tree may return same j via
                        // multiple images for big particles)
                        bool dup = false;
                        for (std::uint32_t kk = 1; kk <= cnt; ++kk) {
                            if (slot_i[kk] == static_cast<std::uint32_t>(j)) {
                                dup = true; break;
                            }
                        }
                        if (dup) continue;
                        if (cnt < cap_i) {
                            ++cnt;
                            slot_i[cnt] = static_cast<std::uint32_t>(j);
                        } else {
                            #pragma omp critical(warn_neighbor_overflow_kd_A)
                            if (!warningIssued) {
                                std::cerr << "Error: particle " << i
                                          << " exceeded MaxNb_i = " << cap_i
                                          << " (kdtree own-write).\n";
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    } else if (refresh[j] == 0) {
                        // Cycle 7 H11: defer cross-write to per-thread
                        // buffer. The serial merge after the parallel
                        // region writes (i) to slot_j without locks; the
                        // sort+unique pass on dirty lists handles dedup.
                        xwrite_buf.emplace_back(
                            static_cast<std::uint32_t>(j),
                            static_cast<std::uint32_t>(i));
                    }
                    // else: j < i and refresh[j] == 1 — j's own walk handles it
                    if (kd_probe_on) tls_kd_write += __rdtsc() - tsc_w0;
                }

                // Advance to next image index combo
                for (std::size_t d = 0; d < Ndim; ++d) {
                    if (++idx[d] < img.counts[d]) break;
                    idx[d] = 0;
                }
            }
        }

        // kdtree Cycle 4: fold thread-locals into atomics on parallel-exit.
        if (kd_probe_on) {
            tls_kd_total = __rdtsc() - tls_kd_tsc_thread_start;
            g_kd_image_cycles.fetch_add(tls_kd_image, std::memory_order_relaxed);
            g_kd_query_cycles.fetch_add(tls_kd_query, std::memory_order_relaxed);
            g_kd_mic_cycles.fetch_add(tls_kd_mic,     std::memory_order_relaxed);
            g_kd_write_cycles.fetch_add(tls_kd_write, std::memory_order_relaxed);
            g_kd_particles_queried.fetch_add(tls_kd_parts, std::memory_order_relaxed);
            g_kd_candidates.fetch_add(tls_kd_cand, std::memory_order_relaxed);
            g_kd_accepted.fetch_add(tls_kd_acc,   std::memory_order_relaxed);
            int kd_tid = omp_get_thread_num();
            if (kd_tid >= 0 && kd_tid < RCP_PROBE_MAX_THREADS) {
                g_kd_thread_cycles[kd_tid].fetch_add(tls_kd_total,
                                                     std::memory_order_relaxed);
            }
        }
    }

    // Cycle 12 H50: sharded parallel cross-write merge tested + reverted.
    // The OMP region setup + per-thread filter-read of ALL buffers
    // cancelled the parallel-write savings — got at-noise +12.5s at
    // N=500K. The serial merge below is well-tuned (~3-4s wall) and
    // simpler. Keep serial.
    for (int t = 0; t < kd_max_threads; ++t) {
        const auto& buf = g_kd_xwrite_buffers[t];
        for (const auto& jp : buf) {
            const std::uint32_t j = jp.first;
            const std::uint32_t i_writer = jp.second;
            std::uint32_t* slot_j = pairs_data + pair_offsets[j];
            const std::uint32_t cap_j = max_neighbors_per_particle[j];
            std::uint32_t cntj = slot_j[0];
            if (cntj < cap_j) {
                slot_j[0] = cntj + 1;
                slot_j[cntj + 1] = i_writer;
                if (!dirty_kd[j]) {
                    kd_dirty_list.push_back(j);
                    dirty_kd[j] = 1;
                }
            } else {
                if (!warningIssued) {
                    std::cerr << "Error: particle " << j
                              << " exceeded MaxNb_j = " << cap_j
                              << " (kdtree cross-write).\n";
                    std::exit(EXIT_FAILURE);
                }
            }
        }
    }

    // Cycle 7 measurement: record query-region wall per regime.
    if (kd_probe_on) {
        const std::uint64_t qregion = __rdtsc() - kd_tsc_qregion_start;
        const int regime = rcp_kd_regime_for_step(g_shadow_current_step);
        g_kd_regime_query_cyc[regime].fetch_add(qregion,
                                                std::memory_order_relaxed);
    }

    // ----- Sort each touched particle's neighbor list -----------------------
    // Cycle 6 H1: skip sort for particles whose lists weren't modified this
    // call. dirty_kd[i] is 1 iff refresh[i] was set (own-write) or a
    // cross-write appended to slot_i.
    // Cycle 7 H11: also dedup via std::unique after sort, since the
    // deferred-buffer merge can introduce duplicates (multi-image returns
    // produce same (j, i) from multiple shifted queries).
    const std::uint64_t kd_tsc_sort_start = kd_probe_on ? __rdtsc() : 0;
    // Cycle 10 H31: iterate the sparse dirty list instead of all N. With
    // mid-regime dirty count ~5K and N=500K, this collapses the parallel-for
    // chunk count from N/static ≈ 125K/thread to ~1K/thread total — cuts the
    // O(N) `if (!dirty_kd[i]) continue` scan that was paying for itself only
    // on the few percent of particles that actually changed.
    const std::size_t kd_D = kd_dirty_list.size();
    #pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t r = 0; r < kd_D; ++r) {
        const std::size_t i = kd_dirty_list[r];
        std::uint32_t* slot = pairs_data + pair_offsets[i];
        std::uint32_t cnt = slot[0];
        if (cnt > 1) {
            std::sort(slot + 1, slot + 1 + cnt);
            auto new_end = std::unique(slot + 1, slot + 1 + cnt);
            slot[0] = static_cast<std::uint32_t>(new_end - (slot + 1));
        }
    }
    if (kd_probe_on) {
        g_kd_sort_cycles.fetch_add(__rdtsc() - kd_tsc_sort_start,
                                   std::memory_order_relaxed);
    }

    // ----- Aggregate observed max/min like sort+walk ------------------------
    // Cycle 12 H44: scan only the sparse dirty list instead of all N. The
    // caller aggregates per-call values monotonically (agg_max only grows,
    // agg_min only shrinks), so a per-call value derived from the touched
    // subset converges to the true bounds across calls. At N=500K with
    // R_dirty ≈ 5K per call, this collapses an O(N) sweep to O(R) — saves
    // ~2.5s wall over a 10K-step bench.
    if (out_max_observed || out_min_observed) {
        std::uint32_t lmax = 0, lmin = std::numeric_limits<std::uint32_t>::max();
        const std::size_t D_count = kd_dirty_list.size();
        for (std::size_t r = 0; r < D_count; ++r) {
            const std::size_t i = kd_dirty_list[r];
            const std::uint32_t c = pairs_data[pair_offsets[i]];
            if (c > lmax) lmax = c;
            if (c < lmin) lmin = c;
        }
        if (out_max_observed) *out_max_observed = lmax;
        if (out_min_observed) *out_min_observed = lmin;
    }

    ++g_pairs_version;
    (void)walls;

    if (kd_probe_on) {
        g_kd_total_cycles.fetch_add(__rdtsc() - kd_tsc_fn_start,
                                    std::memory_order_relaxed);
        g_kd_calls.fetch_add(1, std::memory_order_relaxed);
    }
}

// Dispatch wrapper. The 5 call sites in run_packing_observed call this
// instead of get_pairs_nd_3 directly. Backend behavior:
//   - SortWalk: production path, identical to before scaffolding
//   - KdTree:   exclusive use of kd-tree backend (currently stub == sort+walk)
//   - Shadow:   both run; sort+walk drives forces; kd-tree output compared
static void get_pairs_dispatch(
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
    std::uint32_t* out_max_observed,
    std::uint32_t* out_min_observed)
{
    const PairsBackend backend = rcp_pairs_backend();
    if (backend == PairsBackend::SortWalk) {
        get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                       max_neighbors_per_particle, pair_offsets, pairs_data,
                       out_max_observed, out_min_observed);
        return;
    }
    if (backend == PairsBackend::KdTree) {
        get_pairs_kdtree(N, Ndim, x, D, box, walls, refresh,
                         max_neighbors_per_particle, pair_offsets, pairs_data,
                         out_max_observed, out_min_observed);
        return;
    }
    // Shadow mode: run both, compare, use sort+walk output for forces.
    // CRITICAL: pairs_data carries state across calls — cross-writes in
    // get_pairs_* mutate the EXISTING list for unflagged particles. So
    // shadow_pairs must start with the same state as pairs_data so both
    // backends see the same prior list before refresh-driven mutation.
    const std::uint64_t total = pair_offsets[N];
    static thread_local std::vector<std::uint32_t> shadow_pairs;
    if (shadow_pairs.size() != total) shadow_pairs.resize(total);
    std::memcpy(shadow_pairs.data(), pairs_data,
                total * sizeof(std::uint32_t));
    std::uint32_t shadow_max = 0;
    std::uint32_t shadow_min = std::numeric_limits<std::uint32_t>::max();
    get_pairs_kdtree(N, Ndim, x, D, box, walls, refresh,
                     max_neighbors_per_particle, pair_offsets, shadow_pairs.data(),
                     &shadow_max, &shadow_min);
    get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                   max_neighbors_per_particle, pair_offsets, pairs_data,
                   out_max_observed, out_min_observed);
    compare_pairs_sets(N, pair_offsets, pairs_data, shadow_pairs.data());
}
// ============================================================================
// end kdtree scaffolding
// ============================================================================

void get_forces_nd_3(
    const std::uint32_t* pairs_data,
    const std::uint64_t* pair_offsets,
    const double* x,
    std::size_t N,
    std::size_t Ndim,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    double* F,
    double& U,
    std::vector<std::vector<double>>& min_dist,
    double& max_min_dist,
    double& Lc,
    double& Fmean,
    double mu,
    double& dkappa,
    std::vector<std::size_t>& z)
{
    double K = 1.0;

    // Cycle 7: F is now flat (caller supplies a length-N*Ndim buffer).
    // Parallel zero — same pattern as the per-thread buffers.
    const std::size_t F_size = N * Ndim;
    #pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < F_size; ++k) F[k] = 0.0;
    (void)min_dist;
    if (z.size() == N) {
        std::fill(z.begin(), z.end(), 0);
    } else {
        z.assign(N, 0);
    }
    U = 0.0;
    Lc = 0.0;
    std::size_t count = 1;

    std::vector<double> dx(Ndim);

    max_min_dist = 0;
    bool circle_flag = (walls[0] < 0);

    bool no_walls = true;
    for (std::size_t d = 0; d < Ndim; ++d) {
        if (walls[d] != 0) {
            no_walls = false;
        }
    }

    dkappa = 0;
    Fmean = 0;

    // Cycle 14 layout (restored after Cycle 15b revert): per-thread F_local_flat
    // sized [nthreads * N * Ndim]. F[i] writes go to F_local_flat[tid * stride
    // + i*Ndim + d]; reduce sums across threads into shared F. No atomics.
    //
    // Cycle 15f: pad each per-thread stride up to a 64-byte cache-line
    // boundary so adjacent thread slices don't share a cache line at the
    // seam. Test of the false-sharing hypothesis on T-scaling overhead.
    const int nthreads = omp_get_max_threads();
    constexpr std::size_t CLINE_DBL = 8;       // 64 bytes / sizeof(double)
    constexpr std::size_t CLINE_SZT = 8;       // 64 bytes / sizeof(size_t)
    auto round_up = [](std::size_t v, std::size_t m) {
        return ((v + m - 1) / m) * m;
    };
    const std::size_t F_stride = round_up(N * Ndim, CLINE_DBL);
    const std::size_t z_stride = round_up(N,         CLINE_SZT);
    static std::vector<double>      F_local_flat;
    static std::vector<std::size_t> z_local_flat;
    // Cycle 15f: oversize by 8 doubles / 8 size_t (= 64 bytes) so we can
    // shift the start of slice 0 to a 64-byte aligned address. Combined
    // with cache-line-aligned strides, every per-thread slice starts on a
    // fresh cache line — no boundary sharing between adjacent slices.
    constexpr std::size_t ALIGN_BYTES = 64;
    const std::size_t F_total = static_cast<std::size_t>(nthreads) * F_stride + CLINE_DBL;
    const std::size_t z_total = static_cast<std::size_t>(nthreads) * z_stride + CLINE_SZT;
    g_t_forces_alloc.begin();
    if (F_local_flat.size() != F_total) F_local_flat.resize(F_total);
    if (z_local_flat.size() != z_total) z_local_flat.resize(z_total);
    auto align_up = [](void* p, std::size_t align) -> void* {
        std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(p);
        std::uintptr_t aligned = (addr + align - 1) & ~(static_cast<std::uintptr_t>(align - 1));
        return reinterpret_cast<void*>(aligned);
    };
    double*      F_base = static_cast<double*>(align_up(F_local_flat.data(), ALIGN_BYTES));
    std::size_t* z_base = static_cast<std::size_t*>(align_up(z_local_flat.data(), ALIGN_BYTES));
    #pragma omp parallel num_threads(nthreads)
    {
        const int tid = omp_get_thread_num();
        double* Ft = F_base + tid * F_stride;
        std::size_t* zt = z_base + tid * z_stride;
        std::fill(Ft, Ft + F_stride, 0.0);
        std::fill(zt, zt + z_stride, std::size_t{0});
    }
    g_t_forces_alloc.end();

    // Cycle 14: x_flat shadow refresh eliminated — x is already flat.
    g_t_forces_xflat.begin();
    g_t_forces_xflat.end();

    // Cycle 13: pflat refresh eliminated — pairs are already flat.
    g_t_forces_pflat.begin();
    g_t_forces_pflat.end();

    double U_red = 0.0;
    double Lc_red = 0.0;
    double dkappa_red = 0.0;
    double Fmean_red = 0.0;
    std::size_t count_red = 0;
    double max_min_red = 0.0;

    double*           F_data = F_base;
    std::size_t*      z_data = z_base;
    const double*     x_data       = x;
    const double*     D_data       = D.data();
    const double*     box_data     = box.data();
    // Cycle 9: precompute 1/box per dim so the per-pair MIC wrap becomes
    // floor(delta * inv_box + 0.5) * box instead of floor(delta / box + 0.5)
    // * box. Replaces a division per pair-dim with a multiplication.
    double inv_box[8] = {0.0};
    for (std::size_t d = 0; d < Ndim; ++d) inv_box[d] = 1.0 / box_data[d];

    g_t_forces_loop.begin();
    const bool probe = rcp_cycle_probe_enabled();
    // Cycle 15c MIC-cost probe. Each extra unit injects one extra
    // `floor(delta*inv_box + 0.5) * box` computation per dim, scaled by 1e-16
    // before adding to delta — runtime nonzero so the compiler can't elide,
    // but small enough to be ~1 ulp noise on a delta ~ O(1). Slope of
    // cyc/cand vs extra_mic = per-MIC cost.
    const int extra_mic = rcp_extra_mic_count();
    const double mic_scale = 1e-16;
    #pragma omp parallel default(none) \
        shared(x_data, pairs_data, pair_offsets, D_data, box_data, inv_box, \
               F_data, z_data, N, Ndim, K, F_stride, z_stride, \
               g_floop_cycles, g_floop_candidates, g_floop_overlaps, \
               g_floop_test_cycles, g_floop_overlap_cycles, probe, \
               extra_mic, mic_scale) \
        reduction(+:U_red, Lc_red, dkappa_red, Fmean_red, count_red) \
        reduction(max:max_min_red)
    {
        const int tid = omp_get_thread_num();
        double*      F_t = F_data + tid * F_stride;
        std::size_t* z_t = z_data + tid * z_stride;
        std::uint64_t tls_cycles = 0;
        std::uint64_t tls_candidates = 0;
        std::uint64_t tls_overlaps = 0;
        // Cycle 15c: split per-candidate cycles into "test path" (always
        // executed: load xj, MIC wrap, d2, compare) and "overlap path"
        // (sqrt + F/z buffer writes). Lets us see which side dominates
        // and whether the 38 cyc/cand surplus is in test or overlap.
        std::uint64_t tls_cyc_test = 0;
        std::uint64_t tls_cyc_overlap = 0;
        std::uint64_t tsc_start = probe ? __rdtsc() : 0;

        #pragma omp for schedule(static)
        for (std::size_t i = 0; i < N; ++i) {
            const std::uint32_t* pi = pairs_data + pair_offsets[i];
            const double* xi = x_data + i * Ndim;
            const double Di = D_data[i];
            std::uint32_t numNbr = pi[0];
            // Cycle 6: pairs[i] now stores only neighbors j > i (asymmetric
            // list). The `if (j <= i) continue;` filter is no longer needed —
            // every entry contributes work.
            for (std::uint32_t jdx = 1; jdx <= numNbr; ++jdx) {
                std::size_t j = pi[jdx];

                double r_ij = 0.5 * (Di + D_data[j]);
                const double* xj = x_data + j * Ndim;
                double d2 = 0.0;
                double dx_local[8];

                for (std::size_t d = 0; d < Ndim; ++d) {
                    double delta = xj[d] - xi[d];
                    // Cycle 16d: nearbyint compiles to single vroundsd
                    // (saves the +0.5 op and any sign-handling that floor()
                    // might need). Same MIC math, slightly cheaper.
                    delta -= std::nearbyint(delta * inv_box[d]) * box_data[d];
                    for (int k = 0; k < extra_mic; ++k) {
                        delta += mic_scale *
                                 std::nearbyint(delta * inv_box[d]) * box_data[d];
                    }
                    dx_local[d] = delta;
                    d2 += delta * delta;
                }

                // Cycle 14 H53: count every pair iteration and every applied
                // force. Ratio = rejection rate of the shell-margin pairs
                // (the ones in list but with current distance > sum-of-radii).
                ++tls_candidates;

                const double r_ij_sq = r_ij * r_ij;
                if (d2 < r_ij_sq) {
                    ++tls_overlaps;
                    double dist = std::sqrt(d2);
                    z_t[i]++;
                    z_t[j]++;
                    double F_mag = -K * (r_ij / dist - 1.0);
                    double overlap = 1.0 - dist / r_ij;
                    U_red += overlap;
                    ++count_red;
                    Lc_red += r_ij;
                    if (overlap > max_min_red) { max_min_red = overlap; }

                    for (std::size_t d = 0; d < Ndim; ++d) {
                        double fcomp = F_mag * dx_local[d];
                        F_t[i * Ndim + d] += fcomp;
                        F_t[j * Ndim + d] -= fcomp;
                    }

                    dkappa_red = dkappa_red + K * r_ij * (dist - r_ij);
                    Fmean_red = Fmean_red + F_mag;
                }
            }
        }

        // Cycle 14 H53: always accumulate candidates/overlaps so we can
        // compute force.loop rejection rate without RCP_CYCLE_PROBE.
        g_floop_candidates.fetch_add(tls_candidates, std::memory_order_relaxed);
        g_floop_overlaps.fetch_add(tls_overlaps, std::memory_order_relaxed);
        if (probe) {
            tls_cycles = __rdtsc() - tsc_start;
            g_floop_cycles.fetch_add(tls_cycles, std::memory_order_relaxed);
            g_floop_test_cycles.fetch_add(tls_cyc_test, std::memory_order_relaxed);
            g_floop_overlap_cycles.fetch_add(tls_cyc_overlap, std::memory_order_relaxed);
        }
    }
    g_t_forces_loop.end();

    g_t_forces_reduce.begin();
    // Cycle 14: parallel-over-i reduce. Each i sums across all threads. F[i]
    // is written by exactly one thread (the one whose chunk owns i). Order
    // of summation across t = 0..nthreads-1 is fixed for determinism.
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        double Fi_local[8] = {0};
        std::size_t zi = 0;
        for (int tt = 0; tt < nthreads; ++tt) {
            const double*      F_t = F_data + tt * F_stride + i * Ndim;
            const std::size_t* z_t = z_data + tt * z_stride;
            for (std::size_t d = 0; d < Ndim; ++d) Fi_local[d] += F_t[d];
            zi += z_t[i];
        }
        const std::size_t base = i * Ndim;
        for (std::size_t d = 0; d < Ndim; ++d) F[base + d] += Fi_local[d];
        z[i] += zi;
    }
    g_t_forces_reduce.end();

    U = U_red;
    Lc = Lc_red;
    dkappa = dkappa_red;
    Fmean = Fmean_red;
    count = 1 + count_red;
    max_min_dist = max_min_red;

    if (~no_walls) {

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t d = 0; d < Ndim; ++d) {
                if (walls[d] == 1)
                {
                    double r_ij = D[i] / 2;
                    for (std::size_t wall = 0; wall < 2; ++wall) {
                        dx[d] = box[d] * wall - x[i * Ndim + d];
                        double dist = std::abs(dx[d]);
                        if (dist < r_ij) {
                            double F_mag = -2.0 * K * (r_ij / dist - 1.0);
                            U += 2.0 * (1.0 - dist / r_ij);
                            ++count;
                            double fcomp = F_mag * dx[d];
                            F[i * Ndim + d] += fcomp;
                            dkappa = dkappa + 2 * K * r_ij * (dist - r_ij);
                            Fmean = Fmean + F_mag;
                        }
                    }
                }
            }

            if (circle_flag)
            {
                double r_ij = D[i] / 2;
                double d_ij = 0;
                double R = box[0] / 2;
                for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                {
                    d_ij = d_ij + (x[i * Ndim + M] - R) * (x[i * Ndim + M] - R);
                }
                d_ij = std::sqrt(d_ij);
                double delta = R - d_ij;
                if (delta < r_ij)
                {
                    delta = d_ij;
                    d_ij = 0;
                    for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                    {
                        dx[M] = std::abs(R / delta - 1) * (x[i * Ndim + M] - R);
                        d_ij = d_ij + dx[M] * dx[M];
                    }
                    d_ij = sqrt(d_ij);

                    double F_mag = -2.0 * K * std::abs(r_ij / d_ij - 1.0);

                    U = U + (1 - d_ij / r_ij);
                    ++count;

                    for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                    {
                        F[i * Ndim + M] = F[i * Ndim + M] + F_mag * dx[M];
                    }
                    dkappa = dkappa + 2 * K * r_ij * (d_ij - r_ij);
                    Fmean = Fmean + F_mag;
                }
            }
        }
    }

    U = std::pow(U / double(count), 2);
    Lc = Lc / double(count);

    double sumD = 0.0;
    for (double Di : D) {
        sumD += Di;
    }

    dkappa = dkappa + mu * sumD;
    Fmean = -Fmean / count;
}

void adam_update(
    const std::string& method,
    std::size_t N,
    std::size_t Ndim,
    const double* F,
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
    double& v_hat_kappa)
{
    const double inv_b1_corr = 1.0 / (1.0 - std::pow(beta1, double(t)));
    const double inv_b2_corr = 1.0 / (1.0 - std::pow(beta2, double(t)));
    // Cycle 5a: flat layout for m/v/v_update/m_hat/v_hat means one outer
    // pointer chase per particle (F still nested) and 5 contiguous reads/
    // writes through pre-resolved pointers — vs 6 pointer chases per particle
    // in the prior nested layout. The compiler can now vectorize the inner
    // dim loop because the array bases are loop invariants.
    // Cycle 15l: dropped `v_update[base + kk] = vk;` — vestigial write,
    // never read for computation downstream. Saves ~12% of adam_update's
    // memory write traffic (the kernel is bandwidth-bound at T>1).
    #pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t base = k * Ndim;
        for (std::size_t kk = 0; kk < Ndim; ++kk) {
            double Fkd = F[base + kk];
            double mk = beta1 * m[base + kk] - (1.0 - beta1) * Fkd;
            double vk = beta2 * v[base + kk] + (1.0 - beta2) * Fkd * Fkd;
            m[base + kk] = mk;
            v[base + kk] = vk;
            m_hat[base + kk] = mk * inv_b1_corr;
            v_hat[base + kk] = vk * inv_b2_corr;
        }
    }
    (void)v_update;

    m_kappa = beta1 * m_kappa - (1 - beta1) * dkappa;
    v_kappa = beta2 * v_kappa + (1 - beta2) * (dkappa * dkappa);
    v_update_kappa = v_kappa;
    m_hat_kappa = m_kappa / (1.0 - std::pow(beta1, double(t)));
    v_hat_kappa = v_update_kappa / (1.0 - std::pow(beta2, double(t)));

    (void)method;
    (void)dt;
    (void)verlet_drag;
    (void)v_max;
    (void)a;
    (void)v_verlet;
    (void)a_old;
}

PackingResult run_packing(
    const PackingInput& input,
    const PackingConfig& config)
{
    return run_packing_observed(input, config, 0, false, 0, nullptr).first;
}

std::pair<PackingResult, PackingTrace> run_packing_observed(
    const PackingInput& input,
    const PackingConfig& config,
    std::size_t progress_interval,
    bool capture_positions,
    std::size_t trajectory_interval,
    const PackingObserver& observer,
    const PackingRunOptions& options)
{
    std::size_t N = input.positions.size();
    std::vector<double> D = input.diameters;
    std::vector<double> D0;
    std::size_t Ndim = (input.positions.empty() ? 0 : input.positions[0].size());
    // Cycle 14: x and x_last are flat [N*Ndim]. Convert from the caller's
    // nested input once here; convert back to nested at result time.
    std::vector<double> x(N * Ndim, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t d = 0; d < Ndim; ++d)
            x[i * Ndim + d] = input.positions[i][d];
    std::vector<double> x_last = x;

    std::vector<double> box = config.box;
    std::vector<std::int8_t> walls = config.walls;
    std::uint32_t neighborMax = config.neighbor_max;
    std::uint32_t seed = config.seed;
    bool fix_height = config.fix_height;

    if (box.empty()) box.assign(Ndim, 1.0);
    else if (box.size() != Ndim) {
        throw std::runtime_error("Error: --box length mismatch");
    }
    if (walls.empty()) {
        walls.assign(Ndim, 0);
    }
    else if (walls.size() != Ndim) {
        throw std::runtime_error("Error: --walls expects Ndim entries");
    }

    double phi0 = 0.025;
    double phi = phi0;
    double phi_modifier = 1;
    if (walls[0] < 0)
    {
        if (walls[0] == -1) { walls[0] = -2; }
        for (std::size_t k = 1; k < static_cast<std::size_t>(-walls[0]) && k < Ndim; ++k) { walls[k] = -1; box[k] = box[0]; }

        phi_modifier = sphere_volume(1 / 2., -walls[0]);
    }

    if (fix_height)
    {
        box[Ndim - 1] = box[Ndim - 1] * D[0];
    }

    // Cycle 15o: caller can override the dim-specific initial phi.
    if (options.has_initial_phi_target) {
        phi = options.initial_phi_target;
    }
    else if (Ndim == 2) phi = 0.575;
    else if (Ndim == 3) phi = 0.33;
    else if (Ndim == 4) phi = 0.15;
    else if (Ndim == 5) phi = 0.10;
    else if (Ndim > 5) phi = 0.10 / std::pow(2, (Ndim - 4));
    double delta_phi0 = DELTA_PHI0;
    double delta_phi = DELTA_PHI0;
    double dphi = 0.0;
    std::size_t count = 0;
    bool update_flag = false;
    int direction_flag = 0;
    double F_magnitude = 0.0;
    if (!options.fix_diameter) {
        auto scaled_pair = scale_diameters_nd(D, phi * phi_modifier, box, Ndim, fix_height);
        auto D_scaled = std::move(scaled_pair.first);
        double factor = scaled_pair.second;
        if (fix_height) {
            box[Ndim - 1] = box[Ndim - 1] * factor;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i * Ndim + (Ndim - 1)] = x[i * Ndim + (Ndim - 1)] * factor;
            }
        }
        D = D_scaled;
    }
    double coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
    double current_volume = 0.0;
    for (std::size_t i = 0; i < D.size(); ++i) {
        current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
    }

    double box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());
    phi = current_volume / box_volume / phi_modifier;
    phi0 = phi;

    double Dmin = 0;
    auto it = std::min_element(D.begin(), D.end());
    Dmin = *std::min_element(D.begin(), D.end());
    double Dmin_last = Dmin;
    double max_min_dist = 0;

    if (neighborMax == 0) {
        std::size_t idx = (Ndim >= 2 ? std::min<std::size_t>(Ndim - 2, MAX_NEIGHBORS.size() - 1) : 0);
        neighborMax = MAX_NEIGHBORS[idx];
    }

    std::vector<std::uint32_t> refresh(N, 1);

    // [Cycle 19 cleanup, production lock-in]: per-particle neighbor cap.
    //   scaling = Ndim / 3            (isostatic z(Ndim)/z(3D) — fewer in 2D,
    //                                  more in higher D)
    //   omega   = 1 + Ndim/6          (per-particle exponent — surface area
    //                                  scaling: big 2D particles host fewer
    //                                  neighbors than big 3D particles)
    //   base    = user-supplied config.neighbor_max (legacy keyword), or 200
    //   Nb_i    = min(N, scaling × [50 × (D_i/D_min)^omega + base])
    // For 3D: scaling=1, omega=1.5, base=200 — exactly matches the old
    // empirical formula (50·(D/Dmin)^1.5 + 200).
    auto compute_pair_layout = [&](
        std::vector<std::uint32_t>& max_nb_out,
        std::vector<std::uint64_t>& offsets_out)
    {
        double dmin = D[0];
        for (std::size_t i = 1; i < N; ++i) {
            if (D[i] < dmin) dmin = D[i];
        }
        if (dmin <= 0.0) dmin = 1.0;
        const std::uint32_t base = (neighborMax > 0) ? neighborMax : 200;
        const double scaling = static_cast<double>(Ndim) / 3.0;
        const double omega   = 1.0 + static_cast<double>(Ndim) / 6.0;
        max_nb_out.assign(N, 0);
        offsets_out.assign(N + 1, 0);
        std::uint64_t off = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const double ratio = D[i] / dmin;
            const double cap_d = scaling *
                (50.0 * std::pow(ratio, omega) + static_cast<double>(base));
            std::uint32_t cap = static_cast<std::uint32_t>(std::ceil(cap_d));
            if (static_cast<std::uint64_t>(cap) > N) {
                cap = static_cast<std::uint32_t>(N);
            }
            max_nb_out[i] = cap;
            offsets_out[i] = off;
            off += static_cast<std::uint64_t>(cap) + 1;
        }
        offsets_out[N] = off;
        return off;
    };

    std::vector<std::uint32_t> max_neighbors_per_particle;
    std::vector<std::uint64_t> pair_offsets;
    std::uint64_t pairs_total = compute_pair_layout(
        max_neighbors_per_particle, pair_offsets);
    std::vector<std::uint32_t> pairs_data(pairs_total, 0);
    // Cycle 13 diagnostics: track aggregate observed max/min neighbor count.
    std::uint32_t agg_max_neighbors = 0;
    std::uint32_t agg_min_neighbors = std::numeric_limits<std::uint32_t>::max();
    (void)neighborMax;  // legacy variable no longer used

    // Cycle 7: F flat — contiguous [N*Ndim] for cache-friendly access in
    // adam_update and the forces reduction step.
    std::vector<double> F_flat(N * Ndim, 0.0);
    double* F = F_flat.data();
    // min_dist is a legacy parameter to get_forces_nd_3 but is never read
    // downstream (max_min_dist scalar carries the only consumed value). Pass
    // an empty placeholder instead of allocating N × neighborMax doubles.
    // Saves ~3.2 GB per worker at N=20K.
    std::vector<std::vector<double>> min_dist;
    std::vector<std::size_t> z(N, 0);
    double U = 0.0, Lc = 0.0;

    std::string method = METHOD;
    double beta1 = BETA1, beta2 = BETA2, t = 0, R = 1.00005, alpha = ALPHA_MAX, alpha_max = ALPHA_MAX, N_steps = N_STEPS, epsilon = EPSILON;
    double dt = DT;
    double verlet_drag = 2;
    // Cycle 5a: ADAM state arrays (m, v, v_update, m_hat, v_hat) are flat.
    std::vector<double> m_flat(N * Ndim, 0.0);
    std::vector<double> v_flat(N * Ndim, 0.0);
    std::vector<double> v_update_flat(N * Ndim, 0.0);
    std::vector<double> m_hat_flat(N * Ndim, 0.0);
    std::vector<double> v_hat_flat(N * Ndim, 0.0);
    std::vector<std::vector<double>>
        v_max(N, std::vector<double>(Ndim, 0.0)),
        a(N, std::vector<double>(Ndim, 0.0)),
        v_verlet(N, std::vector<double>(Ndim, 0.0)),
        a_old(N, std::vector<double>(Ndim, 0.0));

    double m_kappa = 0, v_kappa = 0, v_update_kappa = 0, m_hat_kappa = 0, v_hat_kappa = 0;

    std::size_t LastPhiUpdate = 0;
    std::size_t LastAlphaUpdate = 0;
    // Cycle 15: history vectors must be sized for the ACTUAL run length.
    // N_STEPS=60000 is a default cap; options.max_steps can exceed it (e.g.
    // 150K convergence runs). Prior code sized to N_STEPS unconditionally
    // and wrote past end starting at step 60001 — surfaced as glibc
    // 'corrupted size vs. prev_size' ~62K bad writes later.
    const std::size_t hist_size = std::max<std::size_t>(
        static_cast<std::size_t>(N_steps),
        (options.max_steps > 0 ? options.max_steps : 0) + 1);
    std::vector<double> U_history(hist_size, 0.0);

    std::vector<double> phi_history(hist_size, 0.0);
    std::vector<double> F_history(hist_size, 0.0);
    // Cycle 16: ADAM convergence study — capture max overlap, mu, alpha,
    // mu_flag per step for offline diagnostic. Same size as phi/F histories.
    std::vector<double> max_overlap_history(hist_size, 0.0);
    std::vector<double> mu_hist(hist_size, 0.0);
    std::vector<double> alpha_hist(hist_size, 0.0);
    std::vector<std::int8_t> mu_flag_hist(hist_size, 0);

    double phi_min = 0.8;
    if (Ndim == 2) { phi_min = 0.76; }
    if (Ndim == 3) { phi_min = 0.45; }
    if (Ndim == 4) { phi_min = 0.15; }
    if (Ndim == 5) { phi_min = 0.1; }
    if (Ndim > 5) { phi_min = 0.1 / std::pow(2, (Ndim - 4)); }

    double dkappa = 1.0;
    double kappa = 1.0;
    // Cycle 19 Step 5: env var override for initial mu (rung 1 starting mu).
    // Default 2.5e-4: empirically optimal for cyc15 N=100K — phi peak lands
    // just at RCP (~0.836 for cyc15) which gives the oscillation clean
    // starting configurations for basin search. Higher mu_1 over-compresses
    // past RCP (creates damage that oscillation can't recover from);
    // lower mu_1 under-compresses (oscillation can't drive enough caging).
    const double initial_mu_default = [](){
        const char* s = std::getenv("RCP_INITIAL_MU");
        if (!s || !s[0]) return 2.5e-4;
        double v = std::atof(s);
        return (v > 0.0) ? v : 2.5e-4;
    }();
    double mu = options.has_mu_override ? options.mu_override : initial_mu_default;
    bool runtime_fix_diameter = options.fix_diameter;
    bool target_phi_locked = false;
    int mu_flag = 1;
    double Fmean = 0.0;
    double mu_change = 0.;

    std::pair<double, std::size_t> phi_max = {0.0, 0};

    D0 = D;

    // Cycle 4a: cache D0-derived quantities so per-step bookkeep can compute
    // Dmin and current_volume in O(1) instead of two O(N) scans (including
    // an expensive pow() call). Since D[i] = D0[i] * kappa each step (when
    // !runtime_fix_diameter), Dmin = D0_min * kappa and current_volume =
    // coeff * S_D0 * kappa^Ndim where S_D0 = sum_i (D0[i]/2)^Ndim.
    auto recompute_d0_cache = [&]() {
        if (D0.empty()) return std::make_pair(0.0, 0.0);
        double dmin = D0[0];
        double s = 0.0;
        for (std::size_t i = 0; i < D0.size(); ++i) {
            if (D0[i] < dmin) dmin = D0[i];
            double r = D0[i] * 0.5;
            double rN = 1.0;
            for (std::size_t k = 0; k < Ndim; ++k) rN *= r;
            s += rN;
        }
        return std::make_pair(dmin, s);
    };
    auto d0_cache = recompute_d0_cache();
    double D0_min_cached = d0_cache.first;
    double S_D0_cached   = d0_cache.second;

    // Cycle 10 identity tracking: particle_origin[k] tells us which original
    // index the particle currently at slot k started at. Updated alongside
    // every reorder. At the end of the run, used to compare against an
    // un-reordered baseline particle-by-particle.
    std::vector<std::size_t> particle_origin(N);
    for (std::size_t k = 0; k < N; ++k) particle_origin[k] = k;

    // Cycle 17 H6 fix: moved accumulated_refresh declaration here (was at
    // line ~3956) so the reorder lambda below can capture it by reference.
    // Without permuting alongside other state, K-batched refresh flags
    // refer to OLD slot indices after a reorder, causing pair-list churn
    // and F-spikes every 2000 steps.
    std::vector<std::uint32_t> accumulated_refresh(N, 0u);
    bool accumulated_has_refresh = false;
    // Cycle 10: periodic spatial reorder of per-particle arrays. We permute
    // every per-particle array so that particles close in space are also at
    // adjacent indices in memory. The inner force loop's `x[j]` reads then
    // tend to fall in the same cache lines (the prefetcher locks on too).
    //
    // The neighbor-list algorithm itself is untouched — we just renumber
    // particles, then trigger a full pair rebuild. After step ~1000 the
    // system is glassy and the layout stays good for many thousands of steps.
    auto reorder_particles_now = [&]() {
        if (N == 0) return;
        // Cycle 8 H17: 3D Morton spatial reorder (was x-axis column sort).
        // Particles close in 3D space end up close in array index, giving
        // the kdtree leaf scan and force loop better cache locality on the
        // indirect x[]/D[] reads.
        auto perm = (Ndim == 3)
            ? sort_indices_by_morton_3d(x.data(), N, Ndim, box.data())
            : sort_indices_by_column_flat(x.data(), N, Ndim, 0);
        // Track the identity through the permutation so a caller can later
        // map "slot k after the run" back to "original index at run start".
        {
            std::vector<std::size_t> tmp(N);
            for (std::size_t k = 0; k < N; ++k) tmp[k] = particle_origin[perm[k]];
            particle_origin = std::move(tmp);
        }
        // Cycle 14: x is flat. Apply perm with a scatter-then-move pattern.
        {
            std::vector<double> tmp(N * Ndim);
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t src = perm[k] * Ndim;
                const std::size_t dst = k * Ndim;
                for (std::size_t d = 0; d < Ndim; ++d) tmp[dst + d] = x[src + d];
            }
            x = std::move(tmp);
        }
        // x_last gets overwritten by the trailing `x_last = x;` — skip permute.
        // 1D arrays.
        {
            std::vector<double> tmp(N);
            for (std::size_t k = 0; k < N; ++k) tmp[k] = D[perm[k]];
            D = std::move(tmp);
        }
        {
            std::vector<double> tmp(N);
            for (std::size_t k = 0; k < N; ++k) tmp[k] = D0[perm[k]];
            D0 = std::move(tmp);
        }
        // Flat N*Ndim arrays (ADAM state).
        auto apply_perm_flat = [&](std::vector<double>& arr) {
            std::vector<double> tmp(N * Ndim);
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t src = perm[k] * Ndim;
                const std::size_t dst = k * Ndim;
                for (std::size_t d = 0; d < Ndim; ++d) tmp[dst + d] = arr[src + d];
            }
            arr = std::move(tmp);
        };
        apply_perm_flat(m_flat);
        apply_perm_flat(v_flat);
        apply_perm_flat(v_update_flat);
        apply_perm_flat(m_hat_flat);
        apply_perm_flat(v_hat_flat);
        // Cycle 17 H6 fix: reorder triggers full pair-list rebuild via
        // `refresh = [1,1,...]` at line below, so K-batched accumulated
        // refresh flags are no longer relevant — zero them. Original code
        // left them pointing at old slot indices, causing pair-list churn
        // and F-spike on the FIRST post-reorder get_pairs call.
        std::fill(accumulated_refresh.begin(), accumulated_refresh.end(), 0u);
        accumulated_has_refresh = false;

        // Cycle 13: D has been permuted, so per-particle MaxNb caps and
        // offsets shift too. Recompute the layout and resize pairs_data.
        std::uint64_t new_total = compute_pair_layout(
            max_neighbors_per_particle, pair_offsets);
        if (pairs_data.size() != new_total) pairs_data.assign(new_total, 0);
        else std::fill(pairs_data.begin(), pairs_data.end(), 0u);

        // Rebuild pair list fresh under the new numbering.
        std::fill(refresh.begin(), refresh.end(), 1u);
        g_t_pairs.begin();
        {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min);
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
        g_t_pairs.end();
        x_last = x;
    };

    alpha = 0.005;
    alpha_max = 0.005;
    if (Ndim == 2) { alpha = 0.005; alpha_max = 0.005; }
    if (Ndim == 3) {
        alpha = 0.0045; alpha_max = 0.0045;
        // Cycle 17 H18: env var to scale initial alpha (tests faster
        // rung 1 compaction).
        const char* alpha_scale = std::getenv("RCP_ALPHA_SCALE");
        if (alpha_scale && alpha_scale[0]) {
            double s = std::atof(alpha_scale);
            if (s > 0.0 && s < 100.0) { alpha *= s; alpha_max *= s; }
        }
    }
    if (Ndim == 4) { alpha = 0.0035; alpha_max = 0.0035; }
    if (Ndim == 5) { alpha = 0.0025; alpha_max = 0.0025; }
    if (Ndim > 5) { alpha = 0.0025; alpha_max = 0.0025; }

    // ===== profile timer reset (instrumentation) =====
    g_t_pairs.reset(); g_t_forces.reset(); g_t_adam.reset();
    g_t_position.reset(); g_t_refresh.reset(); g_t_bookkeep.reset();
    g_t_forces_alloc.reset(); g_t_forces_xflat.reset();
    g_t_forces_pflat.reset(); g_t_forces_loop.reset();
    g_t_forces_reduce.reset();
    g_t_bk_dupdate.reset(); g_t_bk_meanforce.reset();
    g_t_bk_alpha.reset(); g_t_bk_phistory.reset();
    g_t_bk_misc.reset();
    rcp_cycle_probe_reset();
    // Reset version counters so the pflat-skip optimization doesn't see
    // stale state from a prior run in the same process.
    g_pairs_version = 0;
    g_pflat_version = 0;
    auto profile_wall_start = std::chrono::steady_clock::now();

    // Cycle 10: initial spatial reorder so the very first force loop sees
    // a cache-friendly layout. Skipped if REORDER_INTERVAL=0 (baseline).
    if (options.reorder_interval > 0) {
        reorder_particles_now();
    } else {
        g_t_pairs.begin();
        {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min);
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
        g_t_pairs.end();
    }

    g_t_forces.begin();
    get_forces_nd_3(pairs_data.data(), pair_offsets.data(), x.data(), N, Ndim, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);
    g_t_forces.end();

    PackingTrace trace;
    auto record_observation = [&](std::size_t step_value, bool capture_sample_positions) {
        PackingObservation observation;
        observation.step = step_value;
        observation.phi = phi / phi_modifier;
        observation.force_magnitude = F_magnitude;
        observation.energy = U;
        observation.max_min_dist = max_min_dist;

        if (capture_sample_positions) {
            trace.steps.push_back(step_value);
            trace.diameters.push_back(D);
            trace.phi.push_back(observation.phi);
            trace.force.push_back(observation.force_magnitude);
            trace.energy.push_back(observation.energy);
            trace.max_min_dist.push_back(observation.max_min_dist);
            // Cycle 14: x is flat; convert to nested for the trace API.
            {
                std::vector<std::vector<double>> snap(N, std::vector<double>(Ndim, 0.0));
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t d = 0; d < Ndim; ++d)
                        snap[i][d] = x[i * Ndim + d];
                trace.positions.push_back(std::move(snap));
            }
        }

        if (observer) {
            observer(observation, capture_sample_positions ? &trace.positions.back() : nullptr);
        }
    };

    const bool capture_initial =
        capture_positions && trajectory_interval > 0;
    const bool report_initial =
        observer && progress_interval > 0;
    if (capture_initial || report_initial) {
        record_observation(0, capture_initial);
    }

    std::size_t steps = 0;
    const std::size_t max_steps = options.max_steps > 0 ? options.max_steps : static_cast<std::size_t>(N_steps);
    // Cycle 10: how often to re-reorder. 0 disables (baseline / debugging).
    const std::size_t REORDER_INTERVAL = options.reorder_interval;
    // Cycle 13 b6.A Step 1 (instrumentation) — served its purpose:
    // characterized per-step displacement vs Dmin/4, surfaced the wrap-coord
    // bug in refresh_check (logged box-scale displacements when actual
    // motion was tiny). Kept opt-in via RCP_B6A_PROBE for future audits.
    const bool b6a_probe_on = [](){
        const char* s = std::getenv("RCP_B6A_PROBE");
        return s && s[0] && s[0] != '0';
    }();
    std::vector<double> b6a_x_prev_step;
    if (b6a_probe_on) b6a_x_prev_step = x;
    // Cycle 13 b6.A Phase 1: K-batched get_pairs_dispatch. Refresh flags
    // accumulate across K steps; the expensive range_query work runs only
    // every K steps. Tree correctness is preserved because the patch path
    // at the K-step boundary processes accumulated refresh (motion + wraps)
    // and updates bboxes for all touched leaves. Forces between K-step
    // calls iterate the slightly-stale neighbor list — shell `t` margin
    // absorbs the K-1 step drift.
    //
    // K defaults to 1 (current behavior). Set RCP_K_BATCH=2/3/... to enable.
    const int K_BATCH = [](){
        const char* s = std::getenv("RCP_K_BATCH");
        if (!s || !s[0]) return 1;
        int v = std::atoi(s);
        return (v >= 1 && v <= 16) ? v : 1;
    }();
    // Cycle 18: per-particle troublemaker instrumentation. Tracks aggregate
    // stats over a sliding window (default 5000 steps). When the window
    // fills, accumulators reset and a fresh window starts. At end of run,
    // result holds stats from the LAST window — which is what we care about
    // for identifying particles stuck in oscillation at the end.
    // Cycle 19 Step 5 production lock-in: per-particle troublemaker
    // tracking is research instrumentation only (never feeds back into
    // dynamics). Default OFF; can be re-enabled via RCP_PP_WINDOW > 0
    // for diagnostic runs.
    const std::size_t pp_window_len = [](){
        const char* s = std::getenv("RCP_PP_WINDOW");
        if (!s || !s[0]) return std::size_t(0);
        long v = std::atol(s);
        return (v > 0 && v < 1000000) ? std::size_t(v) : std::size_t(0);
    }();
    const bool pp_enabled = pp_window_len > 0;
    std::vector<double> pp_F_sq_sum(pp_enabled ? N : 0, 0.0);
    std::vector<std::uint64_t> pp_F_sign_flips(pp_enabled ? N : 0, 0);
    std::vector<double> pp_pos_min_run(pp_enabled ? N * Ndim : 0, 0.0);
    std::vector<double> pp_pos_max_run(pp_enabled ? N * Ndim : 0, 0.0);
    std::vector<std::uint64_t> pp_contact_sum(pp_enabled ? N : 0, 0);
    std::vector<std::int8_t> pp_F_sign_prev(pp_enabled ? N * Ndim : 0, 0);
    std::size_t pp_window_steps = 0;
    // Cycle 18 Method B: per-DOF sign-flip momentum reset. If F[k] * m[k] > 0
    // (sign flip vs accumulated EMA direction) AND |F[k]| > threshold, set
    // m[k]=0 BEFORE adam_update. Kills momentum on overshooting DOFs.
    // Enabled via env. Threshold gates against numerical noise.
    const bool sfr_enabled = [](){
        const char* s = std::getenv("RCP_SIGN_FLIP_RESET");
        return s && s[0] && s[0] != '0';
    }();
    const double sfr_F_thresh = [](){
        const char* s = std::getenv("RCP_SIGN_FLIP_F_THRESH");
        if (!s || !s[0]) return 1e-10;  // very low default — let it fire
        double v = std::atof(s);
        return (v > 0.0 && v < 1.0) ? v : 1e-10;
    }();
    std::vector<std::uint64_t> pp_signflip_resets(pp_enabled ? N : 0, 0);
    // Cycle 18 Phase A: F-spike enumeration. Maintain EMA of F (warm-up
    // 200 steps), log any step where F > spike_ratio * EMA. Per spike, log
    // 8 doubles: step, F_ratio, F, maxOv, mu_flag, alpha, step%reorder, step%500.
    // Env-var: RCP_SPIKE_RATIO (default 5.0) sets trigger; 0 disables.
    // Cycle 19 Step 5 production lock-in: F-spike enumeration is research
    // instrumentation only. Default OFF (ratio 0); re-enable via
    // RCP_SPIKE_RATIO > 1 for diagnostic runs.
    const double spike_ratio = [](){
        const char* s = std::getenv("RCP_SPIKE_RATIO");
        if (!s || !s[0]) return 0.0;
        double v = std::atof(s);
        return (v > 1.0 && v < 1000.0) ? v : 0.0;
    }();
    const bool spike_log_enabled = spike_ratio > 0.0;
    std::vector<double> spike_log_buf;
    double F_ema = -1.0;  // initialized on first sample
    const double ema_alpha = 0.005;  // ~200-step memory
    // Cycle 18 Phase A2: reorder event log (always-on at rung -1).
    std::vector<double> reorder_log_buf;
    // Cycle 18 Step 7 Method D: per-particle alpha multiplier for explicit
    // damping of identified troublemakers. Loaded from a file (one line per
    // entry: "particle_id multiplier"). Default mult = 1.0 for all. Only
    // active at rung -1.
    std::vector<double> pp_alpha_mult(N, 1.0);
    bool pp_alpha_mult_loaded = false;
    {
        const char* path = std::getenv("RCP_PP_ALPHA_MULT_FILE");
        if (path && path[0]) {
            std::FILE* fp = std::fopen(path, "r");
            if (fp) {
                std::size_t idx;
                double mult;
                int n_loaded = 0;
                while (std::fscanf(fp, "%zu %lf", &idx, &mult) == 2) {
                    if (idx < N && mult > 0.0 && mult <= 1.0) {
                        pp_alpha_mult[idx] = mult;
                        ++n_loaded;
                    }
                }
                std::fclose(fp);
                pp_alpha_mult_loaded = (n_loaded > 0);
                std::cerr << "[pp_alpha_mult] loaded " << n_loaded
                          << " entries from " << path << "\n";
            } else {
                std::cerr << "[pp_alpha_mult] could not open " << path << "\n";
            }
        }
    }
    // Cycle 19 Step 2: adaptive feedback mu decay. On rung transition, set
    // mu_decay_target = mu/10. Each step, if the system has settled
    // (phi_range over last 250 steps < FLAT_THRESH), tick mu by FACTOR.
    // Self-paces with the system's natural response time: fast systems
    // tick fast, stiff systems tick slowly. No hardcoded period.
    bool mu_decay_active = false;
    double mu_decay_target = 0.0;
    // Cooldown between ticks: enforces a minimum settle window after each
    // tick (else two consecutive ticks on the same flat window would burn
    // through several ticks while the system thinks it's "still flat").
    std::size_t mu_decay_next_tick_step = 0;
    const double mu_decay_factor = [](){
        const char* s = std::getenv("RCP_MU_DECAY_FACTOR");
        if (!s || !s[0]) return 0.85;
        double v = std::atof(s);
        return (v > 0.1 && v < 1.0) ? v : 0.85;
    }();
    const double mu_decay_flat_thresh = [](){
        const char* s = std::getenv("RCP_MU_DECAY_FLAT");
        if (!s || !s[0]) return 5e-3;
        double v = std::atof(s);
        return (v > 0.0) ? v : 5e-3;
    }();
    const std::size_t mu_decay_cooldown = [](){
        const char* s = std::getenv("RCP_MU_DECAY_COOLDOWN");
        if (!s || !s[0]) return std::size_t(250);
        long v = std::atol(s);
        return (v > 0 && v < 100000) ? std::size_t(v) : std::size_t(250);
    }();
    // For the legacy 0->-1 transition trigger which previously needed
    // (step > mu_change + mu_decay_period + 250): keep a coarse upper bound
    // so the transition gate logic still has a "minimum dwell" reference.
    const std::size_t mu_decay_min_dwell = 250;

    // Cycle 19 Step 4: adaptive UP/DOWN search at rung -1, replacing cosine
    // modulation. After 0->-1 decay completes, do N cycles of (UP to
    // base*swing, DOWN to base) using the same tick-on-phi-flat mechanism,
    // then FINAL_DESCENT below base (no target) gated on F + maxOv reaching
    // convergence thresholds. Each tick fires when phi has settled to the
    // current mu level. Same machinery as inter-rung adaptive decay.
    // [Cycle 19 cleanup]: SEARCH algorithm (Step 4 — adaptive UP/DOWN ticks)
    // removed. Superseded by OSC (Step 5) which is the production algorithm.
    // Cycle 19 Step 5: amplitude-decaying oscillation algorithm. Unified
    // replacement for separate 0->-1 adaptive decay + cosine modulation.
    // When 1->0 trigger fires: capture t_basin (= step at trigger), set
    // period = 0.6 * t_basin. Oscillate mu around mu_center, where
    // mu_center exponentially decays from mu_high (= rung 0 base) to
    // mu_low (= rung -1 base) over N cycles. After cycles complete:
    // FINAL_DESCENT below mu_low until F + maxOv gate fires.
    enum class OscPhase {
        INACTIVE,        // before 1->0 trigger or osc disabled
        ACTIVE,          // oscillating around decaying mu_center
        FINAL_DESCENT,   // ticking mu down past mu_low
        DONE             // convergence met
    };
    OscPhase osc_phase = OscPhase::INACTIVE;
    double osc_mu_0 = 0.0;            // rung 0 floor (mu when 1->0 decay completes)
    double osc_mu_minus1 = 0.0;       // rung -1 floor (= mu_0 / 10)
    std::size_t osc_start_step = 0;   // step when oscillation begins
    std::size_t osc_period = 0;       // = osc_period_factor * t_basin
    std::size_t osc_t_relax = 0;      // 1->0 decay duration (measured)
    std::size_t osc_final_start_step = 0;
    std::size_t osc_consec_conv = 0;
    // Cycle 19 Step 5: OSC algorithm is the production default.
    // RCP_ADAPTIVE_OSC=0 disables (falls back to legacy decay + cosine
    // path); any other value (or unset) enables the new oscillation.
    const bool osc_enabled = [](){
        const char* s = std::getenv("RCP_ADAPTIVE_OSC");
        if (!s || !s[0]) return true;
        return s[0] != '0';
    }();
    const double osc_swing = [](){
        const char* s = std::getenv("RCP_OSC_SWING");
        // Default 10: same semantics as the working 10x cosine test — peak
        // at mu_ref*10, trough at mu_ref/10. Reference is now time-decaying
        // from rung-0 floor to rung-1 floor (vs. constant at rung-1 in the
        // original test), but the swing ratio is identical.
        if (!s || !s[0]) return 10.0;
        double v = std::atof(s);
        return (v > 1.0 && v < 1000.0) ? v : 10.0;
    }();
    // Cycle 19 Step 5 / production lock-in: speed preset.
    // RCP_SPEED ∈ {immediate, quick, patient, forever} (default: quick).
    // Sets both rung-1 deadline and oscillation cycle count via the
    // formula:    deadline = Y*4000 + X*sqrt(N*Ndim)
    //             cycles   = Y_cycles  (= 1, 2, 3, 4 across presets)
    // Y, X table: immediate (1, 5, 1), quick (2, 9, 2), patient (3, 24, 3),
    //             forever (4, 49, 4).
    // RCP_RUNG1_DEADLINE and RCP_OSC_CYCLES override the preset if set.
    enum class SpeedPreset { Immediate, Quick, Patient, Forever };
    const SpeedPreset speed_preset = [](){
        const char* s = std::getenv("RCP_SPEED");
        if (!s || !s[0]) return SpeedPreset::Quick;
        std::string v(s);
        if (v == "immediate") return SpeedPreset::Immediate;
        if (v == "quick") return SpeedPreset::Quick;
        if (v == "patient") return SpeedPreset::Patient;
        if (v == "forever") return SpeedPreset::Forever;
        return SpeedPreset::Quick;
    }();
    const std::size_t osc_cycles_total = [speed_preset](){
        const char* s = std::getenv("RCP_OSC_CYCLES");
        if (s && s[0]) {
            long v = std::atol(s);
            if (v >= 1 && v < 100) return std::size_t(v);
        }
        switch (speed_preset) {
            case SpeedPreset::Immediate: return std::size_t(1);
            case SpeedPreset::Quick:     return std::size_t(2);
            case SpeedPreset::Patient:   return std::size_t(3);
            case SpeedPreset::Forever:   return std::size_t(4);
        }
        return std::size_t(2);
    }();
    const double osc_period_factor = [](){
        const char* s = std::getenv("RCP_OSC_PERIOD_FACTOR");
        if (!s || !s[0]) return 0.5;
        double v = std::atof(s);
        return (v > 0.0 && v < 100.0) ? v : 0.5;
    }();
    const double osc_f_thresh = [](){
        const char* s = std::getenv("RCP_OSC_F_THRESH");
        if (!s || !s[0]) return 5e-3;
        double v = std::atof(s);
        return (v > 0.0) ? v : 5e-3;
    }();
    const double osc_ov_thresh = [](){
        const char* s = std::getenv("RCP_OSC_OV_THRESH");
        if (!s || !s[0]) return 5e-4;
        double v = std::atof(s);
        return (v > 0.0) ? v : 5e-4;
    }();
    const std::size_t osc_conv_K = [](){
        const char* s = std::getenv("RCP_OSC_CONV_K");
        if (!s || !s[0]) return std::size_t(3);
        long v = std::atol(s);
        return (v > 0 && v < 100) ? std::size_t(v) : std::size_t(3);
    }();
    const std::size_t osc_max_final_steps = [](){
        const char* s = std::getenv("RCP_OSC_MAX_FINAL");
        if (!s || !s[0]) return std::size_t(50000);
        long v = std::atol(s);
        return (v > 0 && v < 10000000) ? std::size_t(v) : std::size_t(50000);
    }();

    // [Cycle 19 cleanup]: old mu_cosine_* (Cycle 18 Step 9 Method E)
    // declarations removed. Superseded by OSC algorithm.
#if 0  // kept-as-reference block: was the old mu_cosine state — disabled
    const bool mu_cosine_enabled = [](){
        const char* s = std::getenv("RCP_MU_COSINE_RUNG1");
        return s && s[0] && s[0] != '0';
    }();
    const double mu_cosine_swing = [](){
        const char* s = std::getenv("RCP_MU_COSINE_SWING");
        if (!s || !s[0]) return 5.0;
        double v = std::atof(s);
        return (v > 1.0 && v < 1000.0) ? v : 5.0;
    }();
    const std::size_t mu_cosine_period = [](){
        const char* s = std::getenv("RCP_MU_COSINE_PERIOD");
        if (!s || !s[0]) return std::size_t(50000);
        long v = std::atol(s);
        return (v > 100 && v < 10000000) ? std::size_t(v) : std::size_t(50000);
    }();
    double mu_base_rung1 = -1.0;
    std::size_t rung1_start_step = 0;
#endif
    // Cycle 18 Phase B: "oracle mode" — force full refresh of every particle
    // every step at rung -1. Guarantees the kdtree has the correct neighbor
    // set, eliminating the pair-drift gap that causes spike events. Slower
    // per-step but correct. Lets us study Method B without confounding from
    // the pair-tracking issue.
    const bool force_refresh_rung1 = [](){
        const char* s = std::getenv("RCP_FORCE_REFRESH_RUNG1");
        return s && s[0] && s[0] != '0';
    }();
    auto pp_reset_window = [&]() {
        if (!pp_enabled) return;
        std::fill(pp_F_sq_sum.begin(), pp_F_sq_sum.end(), 0.0);
        std::fill(pp_F_sign_flips.begin(), pp_F_sign_flips.end(), 0);
        std::fill(pp_contact_sum.begin(), pp_contact_sum.end(), 0);
        std::fill(pp_signflip_resets.begin(), pp_signflip_resets.end(), 0);
        // Reset pos_min/max to current x
        for (std::size_t k = 0; k < N * Ndim; ++k) {
            pp_pos_min_run[k] = x[k];
            pp_pos_max_run[k] = x[k];
        }
        // Reset previous F sign tracking (treat as 0 initially)
        std::fill(pp_F_sign_prev.begin(), pp_F_sign_prev.end(), 0);
        pp_window_steps = 0;
    };
    if (pp_enabled) pp_reset_window();

    // Cycle 17 H10: optional alpha freeze at rung -1. When env var is set,
    // alpha is locked at its rung-1 entry value (the post-/2 transition
    // alpha). All decay/bump logic is bypassed. Hypothesis: with H6 removing
    // reorder perturbation, alpha dynamics aren't needed at rung -1 — they
    // just create the oscillation cycle that drives F spikes.
    const bool freeze_alpha_rung1 = [](){
        const char* s = std::getenv("RCP_FREEZE_ALPHA_RUNG1");
        return s && s[0] && s[0] != '0';
    }();
    double alpha_rung1_frozen = -1.0;  // -1 means not yet captured
    // Cycle 16 Step 4: ADAM moment reset on mu transition. The 0->-1 reset
    // empirically cuts rung-1 step count 42% and wall 60% at N=50K SR=11
    // by clearing stale m_hat/v_hat statistics carried from the previous
    // mu landscape. Default ON; set RCP_NO_ADAM_RESET_ON_TRANSITION=1 to
    // disable for A/B comparison.
    const bool adam_reset_on_transition = [](){
        const char* s = std::getenv("RCP_NO_ADAM_RESET_ON_TRANSITION");
        return !(s && s[0] && s[0] != '0');
    }();
    auto reset_adam_state = [&]() {
        if (!adam_reset_on_transition) return;
        // Reset m only — v is what protects against tiny-F divergence at
        // step 1 post-reset (step = alpha*F/(sqrt(v)+eps); v=0 means eps
        // dominates and step explodes as alpha*F/eps for small-F DOFs).
        // m is what carried stale direction from the old mu landscape —
        // that's the actual harmful piece. Keep v intact.
        std::fill(m_flat.begin(), m_flat.end(), 0.0);
        m_kappa = 0.0;
        m_hat_kappa = 0.0;
        // Do NOT reset t — bias-correction continues from current step
        // count, preserving the v scaling we just kept.
    };
    // Cycle 16: warm-restart from a saved rung-0-endpoint state. Env-var
    // gated to avoid changing the binding ABI. Loads bit-exact ADAM state
    // (m_flat, v_flat, kappa scalars), positions, diameters, schedule
    // scalars (mu, alpha, mu_flag, mu_change), and the full pre-load
    // history prefixes. Loop step counter resumes at start_step + 1.
    std::size_t start_step = 0;
    {
        const char* load_path = std::getenv("RCP_LOAD_STATE");
        if (load_path && load_path[0]) {
            RcpSavedState ls;
            if (!rcp_load_state(load_path, N, Ndim, ls)) {
                throw std::runtime_error("Cycle 16: failed to load state");
            }
            // Overwrite per-particle state.
            x = ls.positions;
            x_last = x;
            D = ls.diameters;
            D0 = ls.D0;
            m_flat = ls.m_flat;
            v_flat = ls.v_flat;
            // d0_cache, D0_min_cached, S_D0_cached depend on D0; recompute.
            d0_cache = recompute_d0_cache();
            D0_min_cached = d0_cache.first;
            S_D0_cached = d0_cache.second;
            // Pair-layout caps (MaxNb_i and pair_offsets) are computed from
            // D's actual range; the initial call at line ~3637 used the
            // pre-load D which has a different scale. Recompute now so the
            // first get_pairs call after load has correct slot allocation.
            {
                std::uint64_t new_total = compute_pair_layout(
                    max_neighbors_per_particle, pair_offsets);
                if (pairs_data.size() != new_total) {
                    pairs_data.assign(new_total, 0);
                } else {
                    std::fill(pairs_data.begin(), pairs_data.end(), 0u);
                }
            }
            // m_hat/v_hat are recomputed from m,v each step → no need to load.
            // Overwrite schedule scalars.
            alpha = ls.alpha;
            mu = ls.mu;
            kappa = ls.kappa;
            m_kappa = ls.m_kappa;
            v_kappa = ls.v_kappa;
            v_update_kappa = ls.v_update_kappa;
            m_hat_kappa = ls.m_hat_kappa;
            v_hat_kappa = ls.v_hat_kappa;
            t = ls.t_adam;
            phi = ls.phi;
            phi_modifier = ls.phi_modifier;
            phi_max = std::make_pair(ls.phi_max_value, ls.phi_max_step);
            mu_change = ls.mu_change;
            LastPhiUpdate = ls.LastPhiUpdate;
            LastAlphaUpdate = ls.LastAlphaUpdate;
            mu_flag = ls.mu_flag;
            for (std::size_t d = 0; d < Ndim && d < 8; ++d) box[d] = ls.box[d];
            // Splice saved history prefixes into the freshly-allocated
            // hist_size-long vectors. Indices >= ls.hist_used remain 0.
            const std::size_t copy_n = std::min<std::size_t>(
                ls.hist_used, phi_history.size());
            std::copy(ls.phi_history.begin(), ls.phi_history.begin() + copy_n,
                      phi_history.begin());
            std::copy(ls.F_history.begin(), ls.F_history.begin() + copy_n,
                      F_history.begin());
            std::copy(ls.U_history.begin(), ls.U_history.begin() + copy_n,
                      U_history.begin());
            start_step = ls.step;
            // Cycle 18: optional alpha override on resume. Lets us A/B test
            // alpha modifications at rung -1 without having to re-do the
            // 0->-1 transition. RCP_RESUME_ALPHA_MULT scales current alpha;
            // RCP_RESUME_ALPHA_ABS sets absolute value (takes precedence).
            {
                const char* abs_s = std::getenv("RCP_RESUME_ALPHA_ABS");
                const char* mult_s = std::getenv("RCP_RESUME_ALPHA_MULT");
                if (abs_s && abs_s[0]) {
                    double v = std::atof(abs_s);
                    if (v > 0.0 && v < 1.0) {
                        std::cerr << "[rcp_load_state] alpha override: "
                                  << alpha << " -> " << v << "\n";
                        alpha = v;
                    }
                } else if (mult_s && mult_s[0]) {
                    double v = std::atof(mult_s);
                    if (v > 0.0 && v < 100.0) {
                        std::cerr << "[rcp_load_state] alpha scaled: "
                                  << alpha << " -> " << (alpha * v) << "\n";
                        alpha = alpha * v;
                    }
                }
            }
            // Optional ADAM moment reset — tests the staleness hypothesis.
            const char* reset = std::getenv("RCP_RESET_ADAM_ON_LOAD");
            if (reset && reset[0] && reset[0] != '0') {
                std::fill(m_flat.begin(), m_flat.end(), 0.0);
                std::fill(v_flat.begin(), v_flat.end(), 0.0);
                m_kappa = 0.0;
                v_kappa = 0.0;
                v_update_kappa = 0.0;
                m_hat_kappa = 0.0;
                v_hat_kappa = 0.0;
                t = 0;  // bias counter restart
                std::cerr << "[rcp_load_state] ADAM moments RESET "
                             "(m=v=0, t=0)\n";
            }
            std::cerr << "[rcp_load_state] resume start_step="
                      << (start_step + 1) << " max_steps=" << max_steps << "\n";
            if (start_step + 1 > max_steps) {
                throw std::runtime_error(
                    "Cycle 16: saved step exceeds max_steps; nothing to do");
            }
        }
    }
    for (std::size_t step = start_step + 1; step <= max_steps; ++step) {
        rcp_shadow_set_step(static_cast<std::uint64_t>(step));
        g_t_bookkeep.begin();
        steps = step;

        t += 1;

        // Cycle 19 Step 2: feedback-paced mu decay. Tick mu down only when
        // the system has settled to the previous mu (phi_range_250 below
        // FLAT_THRESH). Cooldown prevents firing multiple ticks on the
        // same flat window before the system can respond. Per-step cost
        // when active = one phi-range scan over 250 doubles + a compare;
        // the existing transition gates already do this on every step.
        if (mu_decay_active && step >= mu_decay_next_tick_step && step > 250) {
            // Reuse phi_history (already maintained for transition gates).
            // Note: decay block runs at TOP of step, before phi_history[step-1]
            // is written this iteration. Use [step-250, step-1) to avoid the
            // unfilled trailing slot (matches existing transition-gate pattern).
            auto it_lo = phi_history.begin() + (step - 250);
            auto it_hi = phi_history.begin() + (step - 1);
            auto mm = std::minmax_element(it_lo, it_hi);
            const double phi_range = *mm.second - *mm.first;
            if (phi_range < mu_decay_flat_thresh) {
                mu *= mu_decay_factor;
                if (mu <= mu_decay_target) {
                    mu = mu_decay_target;
                    mu_decay_active = false;
                }
                mu_decay_next_tick_step = step + mu_decay_cooldown;
            }
        }

        // Cycle 19 Step 5: amplitude-decaying oscillation.
        // INACTIVE -> ACTIVE detection: when 1->0 adaptive decay completes
        // (mu_flag==0 and mu_decay_active just dropped false), capture
        // mu_0 (= current mu, rung 0 floor), t_relax (= time the 1->0 decay
        // took), set period = osc_period_factor * t_basin, and start the
        // oscillation at step = current step. Note mu_flag is set to -1 so
        // rung -1 logic (Method B, alpha caps, etc.) applies during osc.
        if (osc_enabled && osc_phase == OscPhase::INACTIVE &&
            mu_flag == 0 && !mu_decay_active && mu_change > 0 &&
            step > mu_change) {
            // mu_change was set at the 1->0 trigger step (= t_basin).
            osc_mu_0 = mu;
            osc_mu_minus1 = mu / 10.0;
            osc_t_relax = step - mu_change;
            if (osc_t_relax < 250) osc_t_relax = 250;  // sanity floor
            osc_start_step = step;
            // Period scales with t_basin (= mu_change here, which is the
            // step at which the rung 1->0 trigger fired).
            osc_period = static_cast<std::size_t>(
                osc_period_factor * static_cast<double>(mu_change));
            if (osc_period < 1000) osc_period = 1000;
            osc_phase = OscPhase::ACTIVE;
            mu_flag = -1;
            mu_change = step;
            // alpha cut at the (skipped) 0->-1 transition. Matches the
            // existing legacy default of 2.0 used in the 0->-1 branch.
            alpha /= 2.0;
        }
        // ACTIVE phase: mu = mu_ref(t) * swing^(-sin(2*pi*t/period))
        //   180° phase shift (-sin instead of sin) so cycle begins with
        //   expansion (trough first) and the cycle-2 endpoint is reached
        //   while mu is descending from a peak. mu_ref exponentially
        //   decays from mu_0 (rung 0 floor) to mu_-1 (rung -1 floor) with
        //   timescale t_relax (the measured 1->0 decay duration).
        if (osc_enabled && osc_phase == OscPhase::ACTIVE) {
            const double t_into = static_cast<double>(step - osc_start_step);
            const double cycle_progress =
                t_into / static_cast<double>(osc_period);
            if (cycle_progress >= static_cast<double>(osc_cycles_total)) {
                mu = osc_mu_minus1;
                osc_phase = OscPhase::FINAL_DESCENT;
                osc_final_start_step = step;
                mu_decay_next_tick_step = step + mu_decay_cooldown;
            } else {
                const double phase_in_cycle =
                    cycle_progress - std::floor(cycle_progress);
                // 180° phase shift: -sin(2*pi*phase)
                const double sin_val = -std::sin(
                    2.0 * 3.14159265358979323846 * phase_in_cycle);
                const double mu_ref =
                    (osc_mu_0 - osc_mu_minus1) *
                    std::exp(-t_into / static_cast<double>(osc_t_relax)) +
                    osc_mu_minus1;
                mu = mu_ref * std::pow(osc_swing, sin_val);

                // Cycle 19 Step 5: peak-gate. After the LAST peak (180° phase
                // puts cycle N peak at t = N - 0.25), open the F+maxOv
                // convergence gate. Catches the high-phi state during descent
                // from compression, rather than waiting for FINAL_DESCENT to
                // anneal mu down to mu_-1 (which loses the peak's density).
                const double last_peak_time =
                    static_cast<double>(osc_cycles_total) - 0.25;
                if (cycle_progress >= last_peak_time) {
                    const double F_now = (step >= 1)
                        ? F_history[step - 1] : 0.0;
                    const bool conditions_met =
                        (F_now < osc_f_thresh) &&
                        (max_min_dist < osc_ov_thresh);
                    if (conditions_met) {
                        osc_consec_conv += 1;
                        if (osc_consec_conv >= osc_conv_K) {
                            osc_phase = OscPhase::DONE;
                        }
                    } else {
                        osc_consec_conv = 0;
                    }
                }
            }
        } else if (osc_enabled && osc_phase == OscPhase::FINAL_DESCENT &&
                   step >= mu_decay_next_tick_step && step > 250) {
            auto it_lo = phi_history.begin() + (step - 250);
            auto it_hi = phi_history.begin() + (step - 1);
            auto mm = std::minmax_element(it_lo, it_hi);
            const double phi_range = *mm.second - *mm.first;
            if (phi_range < mu_decay_flat_thresh) {
                const double F_now = (step >= 1) ? F_history[step - 1] : 0.0;
                const bool conditions_met =
                    (F_now < osc_f_thresh) && (max_min_dist < osc_ov_thresh);
                if (conditions_met) {
                    osc_consec_conv += 1;
                    if (osc_consec_conv >= osc_conv_K) {
                        osc_phase = OscPhase::DONE;
                    }
                    // Hold mu — don't tick further.
                } else {
                    osc_consec_conv = 0;
                    mu *= mu_decay_factor;
                }
                if (step - osc_final_start_step > osc_max_final_steps) {
                    osc_phase = OscPhase::DONE;
                }
                mu_decay_next_tick_step = step + mu_decay_cooldown;
            }
        }
        if (osc_enabled && osc_phase == OscPhase::DONE) {
            steps = step;
            break;
        }

        // [Cycle 19 cleanup]: per-step SEARCH and old cosine modulation
        // blocks removed (both superseded by OSC).

        // Cycle 10: periodic spatial reorder. Skips when REORDER_INTERVAL=0.
        if (REORDER_INTERVAL > 0 && step > 1 && (step % REORDER_INTERVAL == 0)) {
            reorder_particles_now();
        }

        g_t_bk_alpha.begin();
        if ((step - LastAlphaUpdate > 500) &&
            (method == "ADAM") &&
            ((phi - phi_history[std::max<std::size_t>(step - 250, 1)]) < 5e-4) &&
            (step > 2500) && step > mu_change + 500) {

            double trend1 = mean(F_history, step - 450, step - 250) / mean(F_history, step - 75, step - 1);
            double trend2 = (phi - phi_history[std::max<std::size_t>(step - 250, 1)]);

            if ((mu_flag < 1 && trend1 < 0.85) || trend2 < -5E-4) {
                LastAlphaUpdate = step;
                alpha /= 1.5;
            }
            // Cycle 17 H2 REVERTED: suppressing bump-up made it worse —
            // alpha races to zero, system stalls. Bump-up and decay are a
            // needed damped feedback loop. Real lever is capping max alpha
            // at rung -1 (see H5 below). Restoring original bump-up.
            else if (mu_flag < 1 && step - LastAlphaUpdate > 1500) {
                trend1 = mean(F_history, step - 200, step - 1) / mean(F_history, step - 3000, step - 2500);
                if (trend1 > 0.5) {
                    LastAlphaUpdate = step;
                    alpha *= 1.15;
                }
            }

            if (mu_flag == 1) {
                alpha = std::max(alpha, alpha_max / 10.0);
            }
        }

        alpha = std::min(alpha * R, alpha_max);

        if (step % 500 == 0) {
            alpha = std::min(alpha_max, 1.1 * alpha);
        }
        // Cycle 17 H10: hard-lock alpha at rung -1 if freeze is on. Bypasses
        // all decay/bump/cap logic for a clean A/B against the oscillating
        // schedule.
        if (freeze_alpha_rung1 && mu_flag == -1 && alpha_rung1_frozen > 0.0) {
            alpha = alpha_rung1_frozen;
        }
        // Cycle 17 H5: at rung -1, cap alpha at a small fraction of alpha_max.
        // Hypothesis: between spikes, alpha grows by 1.1x per 500 steps
        // (~12x over 14K steps), reaching overshoot threshold and triggering
        // the spike. A hard cap prevents alpha from ever reaching overshoot
        // level. Default ratio 1/100; tunable via RCP_RUNG1_ALPHA_CAP_RATIO.
        if (mu_flag == -1) {
            static const double rung1_alpha_cap_ratio = [](){
                const char* s = std::getenv("RCP_RUNG1_ALPHA_CAP_RATIO");
                if (!s || !s[0]) return 0.01;  // alpha_max / 100
                double v = std::atof(s);
                return (v > 0.0 && v <= 1.0) ? v : 0.01;
            }();
            const double cap = alpha_max * rung1_alpha_cap_ratio;
            if (alpha > cap) alpha = cap;
        }
        g_t_bk_alpha.end();

        g_t_bk_dupdate.begin();
        if (!runtime_fix_diameter) {
            for (std::size_t i = 0; i < D0.size(); ++i) {
                D[i] = D0[i] * kappa;
            }
        }
        g_t_bk_dupdate.end();
        // Cycle 4a: O(1) replacements for the two O(N) D-derived scans below.
        // Dmin_lastest = min(D), but D[i] = D0[i] * kappa for all i, so
        //   min(D) = D0_min_cached * kappa.
        // When runtime_fix_diameter is true, D == D0 and kappa is held at 1
        // (see target_phi branch above), so D0_min_cached IS Dmin.
        double Dmin_lastest = D0_min_cached * (runtime_fix_diameter ? 1.0 : kappa);

        coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
        // current_volume = sum_i coeff * (D[i]/2)^Ndim
        //                = coeff * S_D0_cached * kappa^Ndim   (or kappa^Ndim = 1
        //                  when runtime_fix_diameter and kappa is locked at 1).
        double kappa_to_Ndim = 1.0;
        if (!runtime_fix_diameter) {
            for (std::size_t k = 0; k < Ndim; ++k) kappa_to_Ndim *= kappa;
        }
        current_volume = coeff * S_D0_cached * kappa_to_Ndim;

        box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());

        if (fix_height && !runtime_fix_diameter) {
            box[Ndim - 1] = box[Ndim - 1] * kappa;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i * Ndim + (Ndim - 1)] = x[i * Ndim + (Ndim - 1)] * kappa;
            }
        }

        phi = current_volume / box_volume;

        if (options.has_target_phi && !target_phi_locked && !runtime_fix_diameter) {
            const double user_phi = phi / phi_modifier;
            if (user_phi >= options.target_phi) {
                auto target_scaled = scale_diameters_nd(
                    D,
                    options.target_phi * phi_modifier,
                    box,
                    Ndim,
                    fix_height);
                D = std::move(target_scaled.first);
                D0 = D;
                kappa = 1.0;
                mu = 0.0;
                runtime_fix_diameter = true;
                target_phi_locked = true;

                // D0 has been rescaled — refresh the cached quantities.
                d0_cache = recompute_d0_cache();
                D0_min_cached = d0_cache.first;
                S_D0_cached   = d0_cache.second;

                Dmin_lastest = D0_min_cached;  // kappa = 1 here
                coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
                current_volume = coeff * S_D0_cached;  // kappa^Ndim = 1
                phi = current_volume / box_volume;
                // Cycle 13: D was rescaled, so per-particle MaxNb caps shift.
                {
                    std::uint64_t new_total = compute_pair_layout(
                        max_neighbors_per_particle, pair_offsets);
                    if (pairs_data.size() != new_total) pairs_data.assign(new_total, 0);
                    else std::fill(pairs_data.begin(), pairs_data.end(), 0u);
                }
                std::fill(refresh.begin(), refresh.end(), 1u);
                {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min);
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
                x_last = x;
                Dmin = Dmin_lastest;
            }
        }

        phi_history[step - 1] = phi / phi_modifier;

        if (phi > phi_max.first) {
            phi_max = std::make_pair(phi, step);
        }

        if ((step > 500 && mu_flag == 1)) {
            std::vector<double> XX(phi_history.begin() + step - 250, phi_history.begin() + step - 1);
            double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

            // Cycle 17 H12: rung 1->0 trigger denominator is tunable.
            // Default 3 (original). H12 tested 5 at N=50K SR=11: minor speed
            // change but phi dropped 0.002 (under-compacted, can't recover).
            // Kept the env var hook in case higher-N tests show value.
            static const int rung0_trigger_denom = [](){
                const char* s = std::getenv("RCP_RUNG0_TRIGGER_DENOM");
                if (!s || !s[0]) return 3;
                int v = std::atoi(s);
                return (v >= 2 && v <= 20) ? v : 3;
            }();
            // Cycle 19 Step 5: rung-1 deadline. RCP_RUNG1_DEADLINE direct
            // override takes precedence. Otherwise derive from speed preset:
            //   deadline = Y * 4000 + X * sqrt(N * Ndim)
            //   immediate (Y=1, X=5), quick (Y=2, X=9),
            //   patient   (Y=3, X=24), forever (Y=4, X=49).
            const std::size_t rung1_deadline = [&]() -> std::size_t {
                const char* s = std::getenv("RCP_RUNG1_DEADLINE");
                if (s && s[0]) {
                    long v = std::atol(s);
                    if (v > 100 && v < 10000000) return std::size_t(v);
                }
                int Y = 2, X = 9;  // quick default
                switch (speed_preset) {
                    case SpeedPreset::Immediate: Y = 1; X = 5;  break;
                    case SpeedPreset::Quick:     Y = 2; X = 9;  break;
                    case SpeedPreset::Patient:   Y = 3; X = 24; break;
                    case SpeedPreset::Forever:   Y = 4; X = 49; break;
                }
                const double s_term = static_cast<double>(X) *
                    std::sqrt(static_cast<double>(N) *
                              static_cast<double>(Ndim));
                return static_cast<std::size_t>(
                    static_cast<double>(Y) * 4000.0 + s_term);
            }();
            if (delta_XX < 5e-6 || step - phi_max.second > 3500 || (step > static_cast<int>(rung1_deadline))) {
                // Cycle 19 Step 2: arm adaptive feedback mu decay.
                // mu / 10 (= mu_1 / 10): empirically the winning ratio.
                // Keeps mu_1 → mu_0 a decade decay, so t_relax measured
                // here is consistent with the decade decay used in the
                // oscillation reference exp(-t/t_relax).
                mu_decay_target = mu / 10.0;
                mu_decay_active = true;
                mu_decay_next_tick_step = step + mu_decay_cooldown;
                alpha /= 2.0;
                mu_flag = 0;
                mu_change = step;
            }
        }

        // [Cycle 19 cleanup]: Cycle 17 hardcoded alpha step thresholds
        // (firings at step 45K/55K/60K) removed — never reached at typical
        // run lengths with the new rung-1 deadlines (~10-30K).
        // Cycle 19 Step 2: gate 0->-1 on adaptive decay being COMPLETE
        // (mu_decay_active == false) plus phi-flat. Since each tick fires
        // only on phi_flat at FLAT_THRESH=5e-6, decay completion already
        // implies a settled trajectory; the post-completion check confirms
        // it stays settled at the floor mu.
        // Cycle 19 Step 5: when osc_enabled, the 1->0 adaptive decay just
        // completed; osc handles the rung-0 → rung-1 transition via the
        // exponentially-decaying reference (mu_ref). Suppress the legacy
        // 0->-1 trigger and oscillation start block intercepts instead.
        if (!osc_enabled && step > 500 && mu_flag == 0 && !mu_decay_active &&
            step > mu_change + mu_decay_min_dwell) {
            std::vector<double> XX(phi_history.begin() + step - 250, phi_history.begin() + step - 1);
            double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

            // Cycle 19 Step 2: relaxed 0->-1 gate. Post-decay phi at the
            // rung-0 floor creeps at ~2e-7/step, so phi_range_250 settles
            // around 5e-5 — a 5e-6 gate would never fire under adaptive
            // decay. 2e-5 matches the floor-creep scale.
            if (delta_XX < 2e-5) {
                // Cycle 19 Step 2: arm adaptive feedback mu decay for 0->-1.
                mu_decay_target = mu / 10.0;
                mu_decay_active = true;
                mu_decay_next_tick_step = step + mu_decay_cooldown;
                // Cycle 17 H21: configurable alpha cut at 0->-1 transition.
                // Default 2.0 (original). Larger = smaller initial alpha at
                // rung -1 = less polish-phase noise but slower convergence.
                static const double alpha_cut_0_m1 = [](){
                    const char* s = std::getenv("RCP_ALPHA_CUT_RUNG_M1");
                    if (!s || !s[0]) return 2.0;
                    double v = std::atof(s);
                    return (v >= 1.0 && v <= 100.0) ? v : 2.0;
                }();
                alpha /= alpha_cut_0_m1;
                mu_flag = -1;
                mu_change = step;
                // Cycle 17 H10: capture the rung-1-entry alpha for the
                // freeze-mode lock.
                if (freeze_alpha_rung1) {
                    alpha_rung1_frozen = alpha;
                }
                // Cycle 16 Step 4 (REVERTED): reset at 0->-1 looked great
                // in load-test but in fresh-from-scratch the real bottleneck
                // is alpha oscillation at rung -1 (bump-up vs decay logic
                // fighting → F spikes every 6K-18K steps). The load-test
                // saw a between-spike window; fresh-from-scratch hits them.
                // reset_adam_state();  // SUSPENDED pending alpha co-design.
                // Cycle 16: snapshot the rung-0 endpoint for warm-restart
                // testing. State captured AFTER the transition fires, so the
                // resumed run starts cleanly at step+1 with mu_flag=-1 in
                // effect and mu/alpha already reduced. Opt-in via env var.
                const char* save_path = std::getenv("RCP_SAVE_RUNG0_STATE");
                if (save_path && save_path[0]) {
                    RcpSavedState s;
                    s.step = step;
                    s.hist_used = step;  // [0, step) entries are filled
                    s.alpha = alpha;
                    s.mu = mu;
                    s.kappa = kappa;
                    s.m_kappa = m_kappa;
                    s.v_kappa = v_kappa;
                    s.v_update_kappa = v_update_kappa;
                    s.m_hat_kappa = m_hat_kappa;
                    s.v_hat_kappa = v_hat_kappa;
                    s.t_adam = t;
                    s.phi = phi;
                    s.phi_modifier = phi_modifier;
                    s.phi_max_value = phi_max.first;
                    s.phi_max_step = phi_max.second;
                    s.mu_change = mu_change;
                    s.LastPhiUpdate = LastPhiUpdate;
                    s.LastAlphaUpdate = LastAlphaUpdate;
                    s.mu_flag = mu_flag;
                    for (std::size_t d = 0; d < Ndim && d < 8; ++d) s.box[d] = box[d];
                    s.positions = x;
                    s.diameters = D;
                    s.D0 = D0;
                    s.m_flat = m_flat;
                    s.v_flat = v_flat;
                    s.phi_history.assign(phi_history.begin(), phi_history.begin() + step);
                    s.F_history.assign(F_history.begin(), F_history.begin() + step);
                    s.U_history.assign(U_history.begin(), U_history.begin() + step);
                    rcp_save_state(save_path, N, Ndim, s);
                }
            }
        }

        if (Dmin_lastest / Dmin > 1.05) {
            std::fill(refresh.begin(), refresh.end(), 1u);
            g_t_bookkeep.end();
            g_t_pairs.begin();
            {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min);
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
            g_t_pairs.end();
            g_t_bookkeep.begin();
            x_last = x;
            Dmin = Dmin_lastest;
        } else {
            std::fill(refresh.begin(), refresh.end(), 0u);
            g_t_bookkeep.end();
            g_t_refresh.begin();
            const double dmin_quarter = Dmin / 4.0;
            const double dmin_quarter_sq = dmin_quarter * dmin_quarter;
            // Cycle 15m: was serial; now parallel. Each k is independent
            // (writes to its own refresh[k] / x_last[base + d]). update_flag
            // is monotonically false→true; reduction(||:) is safe.
            // H52 note: dd = x - x_last is intentionally unwrapped. PBC
            // wrap events appear here as box-scale motion and fire refresh
            // — which is correct behavior for the kdtree's spatial structure
            // (the wrapped particle's leaf bbox must update). Fixing the
            // trigger in isolation breaks the kdtree spatial invariant.
            // The K-batching design (Phase 1) addresses the per-step cost
            // by batching the expensive range_query work across K steps;
            // the tree-bbox patch happens at the K-step boundary and
            // handles accumulated wraps correctly.
            #pragma omp parallel for reduction(||:update_flag) schedule(static)
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t base = k * Ndim;
                double d2 = 0.0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double dd = x[base + d] - x_last[base + d];
                    d2 += dd * dd;
                }
                if (d2 > dmin_quarter_sq) {
                    refresh[k] = 1;
                    for (std::size_t d = 0; d < Ndim; ++d)
                        x_last[base + d] = x[base + d];
                    update_flag = true;
                }
            }
            g_t_refresh.end();
            g_t_bookkeep.begin();
            // Cycle 18 Phase B: oracle override — at rung -1 with the env
            // var set, force ALL refresh flags to 1 every step. Eliminates
            // the pair-tracking gap diagnosed in Phase A2 (slow drift over
            // thousands of steps leaving pairs uncounted, then reorder
            // catching them and triggering F-spikes).
            if (force_refresh_rung1 && mu_flag == -1) {
                std::fill(refresh.begin(), refresh.end(), 1u);
                update_flag = true;
                // x_last should track current x so future per-step refresh
                // check is consistent (irrelevant since override always fires
                // at rung -1, but keeps invariant tidy).
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < N * Ndim; ++k) {
                    x_last[k] = x[k];
                }
            }
            // Cycle 13 b6.A Phase 1: accumulate per-step refresh into the
            // K-batch buffer. The accumulated buffer is what get_pairs
            // actually consumes when K-step boundary is hit.
            if (update_flag) {
                accumulated_has_refresh = true;
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < N; ++k) {
                    accumulated_refresh[k] |= refresh[k];
                }
            }
            // Trigger: K-step boundary OR Dmin/4 trigger fired AND K==1
            // (preserve existing per-step behavior at K=1).
            const bool b6a_call_pairs =
                accumulated_has_refresh &&
                (K_BATCH == 1 ? update_flag : (step % static_cast<std::size_t>(K_BATCH) == 0));
            if (b6a_call_pairs) {
                g_t_bookkeep.end();
                g_t_pairs.begin();
                {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, accumulated_refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min);
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
                g_t_pairs.end();
                g_t_bookkeep.begin();
                // Cycle 13 b6.A Phase 1: reset the K-batch accumulator now
                // that the full refresh has run.
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < N; ++k) {
                    accumulated_refresh[k] = 0u;
                }
                accumulated_has_refresh = false;
            }
        }
        update_flag = false;

        g_t_bookkeep.end();
        g_t_forces.begin();
        get_forces_nd_3(pairs_data.data(), pair_offsets.data(), x.data(), N, Ndim, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);
        g_t_forces.end();
        g_t_bookkeep.begin();
        g_t_bk_meanforce.begin();
        F_magnitude = compute_mean_force(F, N, Ndim, Lc, z) / Fmean;
        F_history[step - 1] = F_magnitude;
        // Cycle 16: capture the rest of the per-step diagnostic series.
        // mu, alpha, mu_flag may change later in this step (mu/alpha ladder
        // transitions in the bookkeep block below) — capturing BEFORE those
        // transitions records the state under which this step's force was
        // evaluated. That's the right alignment with F_magnitude.
        max_overlap_history[step - 1] = max_min_dist;
        mu_hist[step - 1] = mu;
        alpha_hist[step - 1] = alpha;
        mu_flag_hist[step - 1] = static_cast<std::int8_t>(mu_flag);
        // Cycle 18 Phase A: spike detection. Wait 200 steps for EMA warmup,
        // then trigger on F > ratio * EMA. Log event details.
        if (spike_log_enabled) {
            if (F_ema < 0.0) {
                F_ema = F_magnitude;
            } else if (step <= 200) {
                F_ema = (1.0 - ema_alpha) * F_ema + ema_alpha * F_magnitude;
            } else {
                // step > 200: enable spike detection
                if (F_magnitude > spike_ratio * F_ema && F_ema > 0.0) {
                    spike_log_buf.push_back(static_cast<double>(step));
                    spike_log_buf.push_back(F_magnitude / F_ema);  // ratio
                    spike_log_buf.push_back(F_magnitude);
                    spike_log_buf.push_back(max_min_dist);
                    spike_log_buf.push_back(static_cast<double>(mu_flag));
                    spike_log_buf.push_back(alpha);
                    spike_log_buf.push_back(static_cast<double>(
                        REORDER_INTERVAL > 0 ? step % REORDER_INTERVAL : 0));
                    spike_log_buf.push_back(static_cast<double>(step % 500));
                }
                F_ema = (1.0 - ema_alpha) * F_ema + ema_alpha * F_magnitude;
            }
        }
        // Cycle 18: per-particle troublemaker accumulators. Updated each
        // step on F, x, and the kdtree's per-particle pair counts.
        if (pp_enabled) {
            if (pp_window_steps >= pp_window_len) pp_reset_window();
            const double* Fd = F;
            const double* xd = x.data();
            #pragma omp parallel for schedule(static)
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t base = k * Ndim;
                double Fsq = 0.0;
                std::uint64_t flips = 0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double f = Fd[base + d];
                    Fsq += f * f;
                    // sign-flip count per dim
                    std::int8_t s = (f > 0) ? 1 : (f < 0 ? -1 : 0);
                    if (s != 0 && pp_F_sign_prev[base + d] != 0 &&
                        s != pp_F_sign_prev[base + d]) {
                        ++flips;
                    }
                    pp_F_sign_prev[base + d] = s;
                    // pos min/max
                    double xv = xd[base + d];
                    if (xv < pp_pos_min_run[base + d]) pp_pos_min_run[base + d] = xv;
                    if (xv > pp_pos_max_run[base + d]) pp_pos_max_run[base + d] = xv;
                }
                pp_F_sq_sum[k] += Fsq;
                pp_F_sign_flips[k] += flips;
                // pair count: max_neighbors_per_particle[k] = MaxNb_i (cap),
                // not the actual count. The actual count is stored as the
                // FIRST slot of pairs_data at pair_offsets[k] (Cycle 13).
                pp_contact_sum[k] += pairs_data[pair_offsets[k]];
            }
            ++pp_window_steps;
        }
        g_t_bk_meanforce.end();

        g_t_bookkeep.end();
        // Cycle 18 Method B: per-DOF sign-flip momentum reset. Pre-step
        // before ADAM. m is the negative EMA of F; if current F has the
        // SAME sign as m (F*m > 0), F has reversed relative to recent
        // average — overshoot detected, kill momentum for that DOF only.
        // |F| threshold gates against numerical noise.
        if (sfr_enabled) {
            double* mf = m_flat.data();
            const double* Fd = F;
            #pragma omp parallel for schedule(static)
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t base = k * Ndim;
                std::uint64_t local_resets = 0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double Fkd = Fd[base + d];
                    if (std::abs(Fkd) > sfr_F_thresh &&
                        Fkd * mf[base + d] > 0.0) {
                        mf[base + d] = 0.0;
                        ++local_resets;
                    }
                }
                if (pp_enabled) pp_signflip_resets[k] += local_resets;
            }
        }
        g_t_adam.begin();
        adam_update(method, N, Ndim,
                    F, dkappa, beta1, beta2, t, dt, verlet_drag,
                    m_flat.data(), v_flat.data(), v_max,
                    v_update_flat.data(),
                    m_hat_flat.data(), v_hat_flat.data(),
                    a, v_verlet, a_old,
                    m_kappa, v_kappa, v_update_kappa,
                    m_hat_kappa, v_hat_kappa);
        g_t_adam.end();
        g_t_bookkeep.begin();

        g_t_bk_phistory.begin();
        if (step > 1100) {
            // Cycle 19 Steps 4 & 5: when adaptive search OR oscillation is
            // active, their own F+maxOv-based convergence gates are the
            // authority. Suppress the legacy phi-flatness gate.
            if (!osc_enabled && mu_flag == -1 &&
                F_history[step - 1] < 5.e-3 &&
                max_min_dist > 1e-16 &&
                step > mu_change + 250) {

                std::vector<double> XX(phi_history.begin() + step - 1000, phi_history.begin() + step - 1);
                double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

                // Cycle 17 H22: phi-flatness threshold is env-var tunable.
                // Cycle 19 default relaxed 2.5e-6 -> 1e-5: the tight threshold
                // was for catching slow upward phi creep across an order of
                // magnitude, which isn't the current regime. 1e-5 allows
                // periods down to ~10K in mu modulation while still ensuring
                // phi has physically settled.
                static const double phi_flat_thresh = [](){
                    const char* s = std::getenv("RCP_CONV_PHI_THRESHOLD");
                    if (!s || !s[0]) return 1e-5;
                    double v = std::atof(s);
                    return (v > 0.0 && v < 1.0) ? v : 1e-5;
                }();
                if (delta_XX < phi_flat_thresh) {
                    g_t_bk_phistory.end();
                    break;
                }
            }
        }
        g_t_bk_phistory.end();

        if (!runtime_fix_diameter) {
            kappa = kappa - alpha * m_hat_kappa / (std::sqrt(v_hat_kappa) + epsilon);
        }
        g_t_bookkeep.end();
        g_t_position.begin();
        if (method != "Verlet") {
            // Cycle 14: x is flat. Pure linear pass through x, mh, vh.
            const double* mh = m_hat_flat.data();
            const double* vh = v_hat_flat.data();
            double* xd = x.data();
            const std::size_t total = N * Ndim;
            // Cycle 18 Method D: per-particle alpha multiplier — active only
            // at rung -1 when loaded. Damps explicit troublemaker particles.
            const bool use_pp_mult = pp_alpha_mult_loaded && (mu_flag == -1);
            if (use_pp_mult) {
                const double* ppm = pp_alpha_mult.data();
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < total; ++k) {
                    const std::size_t pi = k / Ndim;
                    xd[k] -= alpha * ppm[pi] * mh[k] / (std::sqrt(vh[k]) + epsilon);
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < total; ++k) {
                    xd[k] -= alpha * mh[k] / (std::sqrt(vh[k]) + epsilon);
                }
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (std::size_t k = 0; k < N; ++k)
                for (std::size_t d = 0; d < Ndim; ++d)
                    x[k * Ndim + d] += v_verlet[k][d] * dt + 0.5 * a_old[k][d] * dt * dt;
        }

        {
            // Cycle 15k: nested (i,d) loop eliminates the per-element
            // `k % Ndim` modulo. Replace std::fmod (libm transcendental,
            // ~30+ cyc) with a single-step branch wrap (~3 cyc): for
            // typical per-step displacement << box, the wrap is at most
            // one box-add or box-subtract. Matches the actual PBC
            // semantics under bounded-step dynamics.
            double* xd = x.data();
            const double* box_data = box.data();
            #pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < N; ++i) {
                const std::size_t base = i * Ndim;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double v = xd[base + d];
                    const double L = box_data[d];
                    if (v >= L)      v -= L;
                    else if (v < 0)  v += L;
                    xd[base + d] = v;
                }
            }
        }
        g_t_position.end();
        g_t_bookkeep.begin();

        U_history[step - 1] = U;

        const bool should_report =
            observer && progress_interval > 0 && (step % progress_interval == 0);
        const bool should_capture =
            capture_positions && trajectory_interval > 0 && (step % trajectory_interval == 0);
        if (should_report || should_capture) {
            record_observation(step, should_capture);
        }

        if (step % 5000 == 0) {
            std::this_thread::yield();
        }

        // ===== periodic profile snapshot (instrumentation) =====
        if (step % 2000 == 0) {
            auto profile_now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(
                profile_now - profile_wall_start).count();
            std::cerr << "\n[profile snapshot step=" << step
                      << " wall=" << elapsed << "s phi=" << (phi / phi_modifier)
                      << " F=" << F_magnitude
                      << " maxOv=" << max_min_dist
                      << " mu=" << mu
                      << " mu_flag=" << mu_flag
                      << " alpha=" << alpha
                      << "]\n";
            print_timing_summary(elapsed, step);
        }

        (void)delta_phi0;
        (void)delta_phi;
        (void)dphi;
        (void)count;
        (void)direction_flag;
        (void)phi0;
        (void)Dmin_last;
        (void)LastPhiUpdate;
        (void)phi_min;
        (void)seed;
        // Cycle 13 b6.A Step 1 probe: per-step max displacement +
        // refresh_count, logged periodically. Only runs when RCP_B6A_PROBE=1.
        // Uses MIC wrap on the displacement so a PBC wrap shows true motion,
        // not the wrapped-coordinate jump.
        if (b6a_probe_on) {
            double max_step_d2 = 0.0;
            std::size_t refresh_count_step = 0;
            const double* xd = x.data();
            const double* xpd = b6a_x_prev_step.data();
            const double* boxd = box.data();
            const std::size_t Nd = N;
            const std::size_t Nm = Ndim;
            #pragma omp parallel for \
                    reduction(max:max_step_d2) reduction(+:refresh_count_step) \
                    schedule(static)
            for (std::size_t k = 0; k < Nd; ++k) {
                double d2 = 0.0;
                const std::size_t base = k * Nm;
                for (std::size_t d = 0; d < Nm; ++d) {
                    double dd = xd[base + d] - xpd[base + d];
                    const double L = boxd[d];
                    // MIC wrap: handle PBC jumps.
                    if (dd > L * 0.5) dd -= L;
                    else if (dd < -L * 0.5) dd += L;
                    d2 += dd * dd;
                }
                if (d2 > max_step_d2) max_step_d2 = d2;
                if (refresh[k] == 1u) ++refresh_count_step;
            }
            const double max_step_disp = std::sqrt(max_step_d2);
            // Log every 50 steps + the first 20.
            if (step <= 20 || step % 50 == 0) {
                std::cerr << "[b6a-probe] step=" << step
                          << " refresh=" << refresh_count_step
                          << " (" << (100.0 * refresh_count_step / N) << "% of N)"
                          << " max_step_disp=" << max_step_disp
                          << " ratio_to_Dmin=" << (max_step_disp / Dmin)
                          << " ratio_to_Dmin/4=" << (max_step_disp / (Dmin/4.0))
                          << " Dmin=" << Dmin
                          << "\n";
            }
            b6a_x_prev_step = x;  // copy for next-step comparison
        }
        g_t_bookkeep.end();
    }

    // ===== profile timing summary (instrumentation) =====
    {
        auto profile_wall_end = std::chrono::steady_clock::now();
        double total_wall = std::chrono::duration<double>(
            profile_wall_end - profile_wall_start).count();
        print_timing_summary(total_wall, steps);
        rcp_cycle_probe_dump();
        rcp_shadow_dump();
    }

    PackingResult result;
    // Cycle 14: convert flat x back to nested for the result API.
    result.positions.assign(N, std::vector<double>(Ndim, 0.0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t d = 0; d < Ndim; ++d)
            result.positions[i][d] = x[i * Ndim + d];
    result.diameters = std::move(D);
    result.particle_origin = std::move(particle_origin);
    result.box = std::move(box);
    result.walls = std::move(walls);
    result.phi_history = std::move(phi_history);
    result.force_history = std::move(F_history);
    // Cycle 16: surface the diagnostic series to Python.
    result.max_overlap_history = std::move(max_overlap_history);
    result.mu_history = std::move(mu_hist);
    result.alpha_history = std::move(alpha_hist);
    result.mu_flag_history = std::move(mu_flag_hist);
    result.energy_history = std::move(U_history);
    // Cycle 18: per-particle troublemaker stats from the LAST window.
    // pos_range_max is computed from pos_min/max arrays.
    if (pp_enabled) {
        std::vector<double> pos_range_max(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) {
            double mx = 0.0;
            for (std::size_t d = 0; d < Ndim; ++d) {
                double r = pp_pos_max_run[k * Ndim + d] - pp_pos_min_run[k * Ndim + d];
                if (r > mx) mx = r;
            }
            pos_range_max[k] = mx;
        }
        result.pp_F_sq_sum = std::move(pp_F_sq_sum);
        result.pp_F_sign_flips = std::move(pp_F_sign_flips);
        result.pp_pos_range_max = std::move(pos_range_max);
        result.pp_contact_sum = std::move(pp_contact_sum);
        result.pp_signflip_resets = std::move(pp_signflip_resets);
        result.pp_window_steps = pp_window_steps;
    }
    // Cycle 18 Phase A: spike log goes through whether pp is enabled or not.
    result.spike_log = std::move(spike_log_buf);
    result.reorder_log = std::move(reorder_log_buf);
    result.steps = steps;
    result.phi = phi;
    result.max_min_dist = max_min_dist;
    result.force_magnitude = F_magnitude;

    const bool final_already_captured =
        capture_positions && trajectory_interval > 0 && steps > 0 && (steps % trajectory_interval == 0);
    const bool final_already_reported =
        observer && progress_interval > 0 && steps > 0 && (steps % progress_interval == 0);
    if ((capture_positions && !final_already_captured) || (observer && !final_already_reported)) {
        trace.steps.push_back(steps);
        trace.diameters.push_back(result.diameters);
        trace.phi.push_back(phi);
        trace.force.push_back(F_magnitude);
        trace.energy.push_back(U);
        trace.max_min_dist.push_back(max_min_dist);
        if (capture_positions) {
            trace.positions.push_back(result.positions);
        }
        if (observer) {
            PackingObservation observation;
            observation.step = steps;
            observation.phi = phi;
            observation.force_magnitude = F_magnitude;
            observation.energy = U;
            observation.max_min_dist = max_min_dist;
            observer(observation, capture_positions ? &trace.positions.back() : nullptr);
        }
    }

    return {std::move(result), std::move(trace)};
}

}  // namespace rcpgenerator
