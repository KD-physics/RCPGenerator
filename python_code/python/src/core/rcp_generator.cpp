#include "rcpgenerator/rcp_generator.hpp"
#include "rcpgenerator/diameter_permute.hpp"
#include "rcpgenerator/ellipsoid.hpp"
#include "rcpgenerator/rod.hpp"
#include "rcpgenerator/curved_wall.hpp"

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

#if defined(__linux__)
#include <sys/mman.h>   // madvise(MADV_DONTNEED): release pairs_data pages on reorder
#endif

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
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__i386__) || defined(__x86_64__)
#include <x86intrin.h>
#endif
#include <cstdlib>

static inline std::uint64_t rcp_read_cycles() {
#if defined(_MSC_VER) || defined(__i386__) || defined(__x86_64__)
    return __rdtsc();
#else
    return 0;
#endif
}

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
// Cycle 21: AVX-512 8-wide batching of the force-loop reject test.
// Default ON where eligible (__AVX512F__ build, Ndim==3, extra_mic==0;
// scalar fallback otherwise). RCP_FORCE_SIMD=0 disables. Lane arithmetic
// mirrors the scalar path operation-for-operation (same FMA contraction
// shape, roundscale == nearbyint under default rounding), and lanes are
// processed in jdx order — VALIDATED bit-exact vs scalar by the
// fixed-state byte A/B (3000 steps: phi/force/maxOv histories and final
// positions/diameters all hash-identical).
static bool rcp_force_simd() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_FORCE_SIMD");
        if (!s || !s[0]) return true;
        return s[0] != '0';
    }();
    return v;
}
// Cycle 21k: packed-D0 mirror for the force loop's reject test. Default
// OFF — measured NET NEGATIVE (+22% wall at 100K): D is cache-resident
// at any practical N (4MB at 500K vs 16MB L3), so the D[j] gather it
// replaces is nearly free, while streaming the cold ~8B/slot mirror
// through the force loop evicts useful cache lines. Kept (RCP_PACKED_D0=1)
// as a documented negative, bit-exact-validated mechanism.
static bool rcp_packed_d0() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_PACKED_D0");
        if (!s || !s[0]) return false;
        return s[0] != '0';
    }();
    return v;
}
// Cycle 21L FP32 STUB (measurement only, default OFF, never ships):
// the 3D SIMD reject test gathers positions from a float shadow and
// converts in-register — same gather element count, half the bytes.
// Isolates the bandwidth half of the FP32-engine estimate before the
// 16-lane build-out. NOT byte-exact (float-rounded coordinates feed the
// test and the dx values), so it is gated and excluded from byte A/Bs.
static bool rcp_fp32_stub() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_FP32_STUB");
        if (!s || !s[0]) return false;
        return s[0] != '0';
    }();
    return v;
}
// Cycle 21M rung-1 mu boost (Ken's design, default OFF): while at rung 1
// (mu_flag == 1), the mu fed to the force/dkappa computation is
//   mu_eff = mu_1 * (1 + (B-1) * exp(-step/tau)),
// i.e. it starts at B x mu_1 and relaxes to mu_1 on timescale tau. The
// drive is front-loaded into the dilute regime (fast quench where rate
// doesn't matter) and hands over to the standard mu_1 well before the
// dense regime where static-high mu_1 was measured to damage the basin
// (over-compression; phi overshoot on relax). The `mu` VARIABLE is never
// modified, so the rung 1->0 transition machinery (mu_decay_target =
// mu/10, oscillation bases) is untouched by construction.
// H12c [163]: reopen mu-kick fraction for variant-B bottleneck (small fraction of the
// rung-1 mu -> slight overlap ~1e-4 that reopens cages a bit). File-scope so it adds no
// static-local guard inside the byte-exact run_packing_observed hot path.
static double rcp_bneck_kick_frac() {
    static const double v = []() {
        const char* s = std::getenv("RCP_BNECK_KICK_FRAC");
        if (!s || !s[0]) return 0.05;
        double x = std::atof(s);
        return (x > 0.0 && x < 1.0) ? x : 0.05;
    }();
    return v;
}
static double rcp_mu1_boost() {
    static const double v = []() {
        const char* s = std::getenv("RCP_MU1_BOOST");
        // Cycle 21M default ON (B=5): statistically validated at 100K
        // SR=35 3D, 3 seeds/arm — phi-neutral-to-positive at the default
        // deadline and REQUIRED for the -2000 deadline cut to be
        // phi-free (cut-only loses 4.3e-4; cut+boost loses 1.0e-4,
        // within run-to-run noise). RCP_MU1_BOOST=1 disables.
        if (!s || !s[0]) return 5.0;
        double x = std::atof(s);
        return (x >= 1.0 && x <= 20.0) ? x : 5.0;
    }();
    return v;
}
static double rcp_mu1_boost_tau() {
    static const double v = []() {
        const char* s = std::getenv("RCP_MU1_BOOST_TAU");
        if (!s || !s[0]) return 350.0;
        double x = std::atof(s);
        return (x >= 10.0 && x <= 10000.0) ? x : 350.0;
    }();
    return v;
}
// Cycle 21k adaptive shell: the run loop widens the pair shell while the
// diameter growth rate is high (early compression), so refresh-trigger
// budgets last longer and the early-phase kdtree query storm collapses.
// RCP_SHELL_SMAX caps the multiplier (1.0 disables and reproduces the
// fixed-shell behavior byte-for-byte); RCP_SHELL_K is the target list
// lifetime in steps under pure-growth budget burn.
// Default 1.0 (disabled) — measured NET NEGATIVE at SMAX=1.5/K=40 even
// with the signed-EMA controller (+3% wall at 100K SR=35): the early
// compression phase starts at phi ~0.55 where rows are already dense, so
// the widened shell's force-loop cost (+6%) outweighs the pairs saving
// (-6% on a component half the size). Mechanism kept (RCP_SHELL_SMAX>1)
// with the per-particle min-clamped budget proof intact.
static double rcp_shell_smax() {
    static const double v = []() {
        const char* s = std::getenv("RCP_SHELL_SMAX");
        if (!s || !s[0]) return 1.0;
        double x = std::atof(s);
        return (x >= 1.0 && x <= 2.0) ? x : 1.0;
    }();
    return v;
}
static double rcp_shell_k() {
    static const double v = []() {
        const char* s = std::getenv("RCP_SHELL_K");
        if (!s || !s[0]) return 40.0;
        double x = std::atof(s);
        return (x >= 1.0 && x <= 1000.0) ? x : 40.0;
    }();
    return v;
}
// Cycle 21c: fuse particle ADAM + position update + PBC wrap into one
// pass (default ON for non-Verlet; RCP_FUSED_ADAM=0 reverts to the split
// kernels). Per-element math is verbatim adam_update + position + wrap;
// m_hat/v_hat live in registers (their only array reader was the position
// pass). Validated by the fixed-state byte A/B.
static bool rcp_fused_adam() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_FUSED_ADAM");
        if (!s || !s[0]) return true;
        return s[0] != '0';
    }();
    return v;
}
static PairsBackend rcp_pairs_backend() {
    // Cycle 20: kdtree promoted to default backend. At N=10K/T=1 it cuts
    // neighbor-management wall time 3.0-3.5x vs sort+walk (2D and 3D), and
    // it passed the missed-pair sentinel on full convergence runs in both
    // dimensions with the growth-aware refresh trigger. sort+walk remains
    // available via RCP_PAIRS_BACKEND=sortwalk (it is also still the
    // trusted reference the sentinel builds against internally).
    static const PairsBackend v = []() {
        const char* s = std::getenv("RCP_PAIRS_BACKEND");
        if (!s || !*s) return PairsBackend::KdTree;
        if (std::strcmp(s, "sortwalk") == 0) return PairsBackend::SortWalk;
        if (std::strcmp(s, "kdtree") == 0) return PairsBackend::KdTree;
        if (std::strcmp(s, "shadow") == 0) return PairsBackend::Shadow;
        return PairsBackend::KdTree;
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
static KernelTimer g_t_swap;      // diameter-SWAP sweeps (zero calls when swap == 0)

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
// Cycle 21d: permutation epoch. Incremented by the spatial reorder so the
// kdtree's D-stat index cache invalidates DETERMINISTICALLY. The previous
// 2-sentinel ratio check (D[0], D[N/2]) can false-negative when the
// deterministic-quantile sampler creates size ties and a permutation maps
// both sentinels onto same-sized particles — caught by the missed-pair
// sentinel as a growing miss window between two reorders (stale Dmax index
// -> query radius too small for big-particle pairs).
static std::uint64_t g_rcp_perm_epoch = 0;
// checkpoint-deterministic mode [log 044/B4]: when a run checkpoints, the
// kd/oct tree must hold NO cross-step state (the patch path's leaf staleness
// is path-dependent and unserializable at reasonable cost), so every pair
// update takes the full-rebuild branch. Statistically equivalent to the
// patched default; bitwise different (documented mode semantics).
static bool g_rcp_ckpt_deterministic = false;
// kdtree D-stat cache, hoisted from get_pairs so checkpoint/resume can
// serialize it: the cache is PATH-DEPENDENT (stale indices survive hits) and
// bitwise resume must reproduce the unbroken run's staleness [log 044/B3].
static std::size_t g_dstat_max_idx = 0, g_dstat_min_idx = 0, g_dstat_med_idx = 0;
static std::size_t g_dstat_N = 0;
static double g_dstat_D0 = 0.0, g_dstat_Dmid = 0.0;
static bool g_dstat_valid = false;
static std::uint64_t g_dstat_epoch = ~0ull;
// Cycle 21k adaptive shell multiplier. Single writer (the run loop, once
// per step before the refresh check); the pair backends read it inside
// the dispatch call that follows on the same thread.
static double g_rcp_shell_scale_dyn = 1.0;

// Cycle 21L diagnostics: current RSS of this process (Linux), GB.
static double rcp_rss_gb() {
    FILE* f = std::fopen("/proc/self/statm", "r");
    if (!f) return 0.0;
    long tot = 0, res = 0;
    if (std::fscanf(f, "%ld %ld", &tot, &res) != 2) res = 0;
    std::fclose(f);
    return res * 4096.0 / 1e9;
}
static void rcp_rss_probe(const char* tag) {
    if (!std::getenv("RCP_DEBUG_LAYOUT")) return;
    std::cerr << "[rss] " << tag << ": " << rcp_rss_gb() << " GB\n";
}
// Neighbor-capacity overflow: build the error message thrown to Python.
// Sites inside OpenMP regions record it under a critical section and the
// enclosing function throws AFTER the region (throwing across an OpenMP
// boundary would call std::terminate and kill the process — which is the
// failure mode this replaces).
static std::string rcp_nb_overflow_msg(const char* site, std::size_t particle,
                                       std::uint32_t cap) {
    std::string m = "rcpgenerator: particle ";
    m += std::to_string(particle);
    m += " exceeded its neighbor capacity (";
    m += std::to_string(cap);
    m += " slots) during pair-list construction (";
    m += site;
    m += ").";
    m += " This usually means the per-particle neighbor allocation is too small for the size ratio in this packing. Recommendation: pass a larger neighbor_max (the per-particle base; e.g. neighbor_max=1500 -- capacity scales as 50*(D_i/D_min)^omega + neighbor_max), or reduce the size ratio.";
    return m;
}
// Per-run/periodic timing diagnostics to stderr. Default OFF for release
// builds; RCP_PROFILE=1 enables (research instrumentation).
static bool rcp_profile_prints() {
    static const bool v = []() {
        const char* s = std::getenv("RCP_PROFILE");
        return s && s[0] && s[0] != '0';
    }();
    return v;
}
static void print_timing_summary(double total_wall, std::size_t steps) {
    if (!rcp_profile_prints()) return;
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
    show("swap_sweep       ", g_t_swap);
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
    std::uint32_t* out_min_observed,
    double* packed_d0,
    const double* D0_data)
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
    t *= g_rcp_shell_scale_dyn;   // Cycle 21k adaptive shell

    const bool probe_sort = rcp_cycle_probe_enabled();
    std::uint64_t sort_tsc0 = probe_sort ? rcp_read_cycles() : 0;
    auto sort_idx = sort_indices_by_column_flat(x, N, Ndim, 0);
    if (probe_sort) {
        g_gpairs_sort_cycles.fetch_add(rcp_read_cycles() - sort_tsc0,
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
    for (std::size_t d = 0; d < Ndim; ++d) gp_inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;  // [235] minimum image only along periodic dimensions: 0 along a wall -> no wrap

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
    std::string overflow_msg;   // set under the overflow criticals below

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
        std::uint64_t tsc_start = probe_gp ? rcp_read_cycles() : 0;
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
                                warningIssued = true;
                                overflow_msg = rcp_nb_overflow_msg(
                                    "sort+walk own-write", i, cap_i);
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
                                    warningIssued = true;
                                    overflow_msg = rcp_nb_overflow_msg(
                                        "sort+walk cross-write", j, cap_j);
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
            tls_cycles = rcp_read_cycles() - tsc_start;
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

    // Cycle 21N: an overflow inside the parallel region was recorded
    // instead of exiting; surface it now as an exception (caught in
    // Python; the worker process survives).
    if (warningIssued) throw std::runtime_error(overflow_msg);

    // Canonicalize: sort each particle's neighbor list.
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        std::uint32_t* slot = pairs_data + pair_offsets[i];
        std::uint32_t cnt = slot[0];
        if (cnt > 1) {
            std::sort(slot + 1, slot + 1 + cnt);
        }
        // Cycle 21k: refill the packed-D0 mirror for this row. Every row
        // this call mutated (refresh rebuilds AND cross-write appends)
        // passes through here, so the mirror can never go stale.
        if (packed_d0) {
            double* pk = packed_d0 + pair_offsets[i];
            for (std::uint32_t k = 1; k <= cnt; ++k)
                pk[k] = D0_data[slot[k]];
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
// Cycle 20: missed-pair sentinel (validation instrument, env-gated).
//
// RCP_SENTINEL=<S> (default 0 = off): every S steps, build a trusted full
// pair set with get_pairs_nd_3 (refresh-all, scratch buffers) and verify
// that every ACTUALLY OVERLAPPING pair — dist < (Di+Dj)/2, i.e. every
// force-generating pair — appears in the live candidate list used for this
// step's forces. This directly tests the hard invariant "no missed true
// interacting pairs" for whatever backend / refresh-trigger / K-batch
// configuration is live, at every phase of the run (early/mid/late).
//
// Each missed pair prints a stderr line with step, indices, and fractional
// overlap. End-of-run summary prints totals. Caveat: per-particle caps can
// truncate the scratch list exactly as they would a live list (warning is
// emitted by get_pairs in that case), so a capped configuration validates
// against the same capped universe it runs with.
// ============================================================================
static long rcp_sentinel_every() {
    static const long v = []() {
        const char* s = std::getenv("RCP_SENTINEL");
        return (s && s[0]) ? std::atol(s) : 0L;
    }();
    return v;
}
static std::uint64_t g_sentinel_checks = 0;
static std::uint64_t g_sentinel_overlap_pairs = 0;
static std::uint64_t g_sentinel_missed = 0;
static double g_sentinel_worst_missed_overlap = 0.0;
static double g_sentinel_seconds = 0.0;

static bool rcp_live_list_has(const std::uint32_t* pairs_data,
                              const std::uint64_t* pair_offsets,
                              std::size_t a, std::uint32_t b) {
    const std::uint32_t* slot = pairs_data + pair_offsets[a];
    const std::uint32_t cnt = slot[0];
    for (std::uint32_t k = 1; k <= cnt; ++k)
        if (slot[k] == b) return true;
    return false;
}

static void rcp_sentinel_check(
    std::size_t N, std::size_t Ndim,
    const double* x, const std::vector<double>& D,
    const std::vector<double>& box, const std::vector<std::int8_t>& walls,
    const std::uint32_t* max_neighbors_per_particle,
    const std::uint64_t* pair_offsets,
    const std::uint32_t* live_pairs_data,
    std::size_t step,
    // Cycle 21i forensics (optional, may be null): trigger-state dump
    // for missed pairs.
    const double* x_last = nullptr,
    const double* D_last = nullptr,
    const std::uint32_t* refresh_flags = nullptr,
    double theta = 0.0)
{
    const auto t0 = std::chrono::steady_clock::now();
    static std::vector<std::uint32_t> scratch_pairs;
    static std::vector<std::uint32_t> refresh_all;
    const std::uint64_t total = pair_offsets[N];
    if (scratch_pairs.size() < total) scratch_pairs.resize(total);
    if (refresh_all.size() != N) refresh_all.assign(N, 1u);

    // Observation-only: restore the pairs-version counter the scratch
    // build bumps (currently unread, but keep the instrument side-effect
    // free by construction).
    const std::uint64_t saved_pairs_version = g_pairs_version;
    std::uint32_t dummy_max = 0, dummy_min = 0;
    get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh_all,
                   max_neighbors_per_particle, pair_offsets,
                   scratch_pairs.data(), &dummy_max, &dummy_min);
    g_pairs_version = saved_pairs_version;

    double inv_box[8] = {0.0};
    for (std::size_t d = 0; d < Ndim; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;  // [235] minimum image only along periodic dimensions: 0 along a wall -> no wrap

    std::uint64_t overlaps = 0, missed = 0;
    double worst = 0.0;
    #pragma omp parallel for schedule(dynamic, 64) \
        reduction(+:overlaps, missed) reduction(max:worst)
    for (std::size_t i = 0; i < N; ++i) {
        const std::uint32_t* slot = scratch_pairs.data() + pair_offsets[i];
        const std::uint32_t cnt = slot[0];
        const double* xi = x + i * Ndim;
        const double Di = D[i];
        for (std::uint32_t k = 1; k <= cnt; ++k) {
            const std::uint32_t j = slot[k];
            const double* xj = x + std::size_t(j) * Ndim;
            const double r_ij = 0.5 * (Di + D[j]);
            double d2 = 0.0;
            for (std::size_t d = 0; d < Ndim; ++d) {
                double delta = xj[d] - xi[d];
                delta -= std::nearbyint(delta * inv_box[d]) * box[d];
                d2 += delta * delta;
            }
            if (d2 >= r_ij * r_ij) continue;   // not force-generating
            ++overlaps;
            if (rcp_live_list_has(live_pairs_data, pair_offsets, i, j) ||
                rcp_live_list_has(live_pairs_data, pair_offsets, j, i))
                continue;
            ++missed;
            const double frac = 1.0 - std::sqrt(d2) / r_ij;
            if (frac > worst) worst = frac;
            #pragma omp critical(rcp_sentinel_print)
            {
                std::cerr << "[sentinel] MISSED PAIR step=" << step
                          << " i=" << i << " j=" << j
                          << " overlap_frac=" << frac
                          << " Di=" << Di << " Dj=" << D[j]
                          << " dist=" << std::sqrt(d2) << "\n";
                // Forensics: live slot counts and caps for both particles.
                const std::uint32_t ci =
                    live_pairs_data[pair_offsets[i]];
                const std::uint32_t cj =
                    live_pairs_data[pair_offsets[j]];
                std::cerr << "[sentinel]   slots: cnt_i=" << ci
                          << "/" << max_neighbors_per_particle[i]
                          << " cnt_j=" << cj
                          << "/" << max_neighbors_per_particle[j] << "\n";
                // Cycle 21i: trigger-state forensics.
                if (x_last && D_last && refresh_flags) {
                    for (const std::size_t k :
                         {static_cast<std::size_t>(i),
                          static_cast<std::size_t>(j)}) {
                        double md2 = 0.0;
                        for (std::size_t d = 0; d < Ndim; ++d) {
                            const double dd =
                                x[k * Ndim + d] - x_last[k * Ndim + d];
                            md2 += dd * dd;
                        }
                        std::cerr << "[sentinel]   k=" << k
                                  << " |dx|=" << std::sqrt(md2)
                                  << " dD/2="
                                  << 0.5 * (D[k] - D_last[k])
                                  << " theta=" << theta
                                  << " refresh=" << refresh_flags[k] << "\n";
                    }
                }
            }
        }
    }
    ++g_sentinel_checks;
    g_sentinel_overlap_pairs += overlaps;
    g_sentinel_missed += missed;
    if (worst > g_sentinel_worst_missed_overlap)
        g_sentinel_worst_missed_overlap = worst;
    g_sentinel_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}

static void rcp_sentinel_dump() {
    if (rcp_sentinel_every() <= 0) return;
    std::cerr << "\n=== Missed-pair sentinel summary ===\n"
              << "  checks:            " << g_sentinel_checks << "\n"
              << "  overlap pairs:     " << g_sentinel_overlap_pairs << "\n"
              << "  MISSED pairs:      " << g_sentinel_missed << "\n"
              << "  worst missed frac: " << g_sentinel_worst_missed_overlap << "\n"
              << "  sentinel cost:     " << g_sentinel_seconds << "s\n"
              << "  VERDICT:           "
              << (g_sentinel_missed == 0 ? "PASS (no missed pairs)"
                                         : "*** FAIL ***") << "\n";
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
            // Cycle 21N FIX: bound the partition to the LIVE N entries.
            // perm is a grow-only buffer; using perm.end() partitioned
            // across stale tail entries from a previous (larger-N) run in
            // the same process, shuffling dead particle indices into the
            // live range — out-of-bounds position reads (occasional
            // worker segfault) and silently wrong neighbor lists. Latent
            // since Cycle 8 H19; only fires when one process packs a
            // smaller system after a larger one.
            __gnu_parallel::nth_element(
                perm.begin(), perm.begin() + mid,
                perm.begin() + static_cast<std::ptrdiff_t>(N),
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
#if defined(__GNUC__) && !defined(__clang__)
    __attribute__((optimize("fast-math", "tree-vectorize")))
#endif
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
#if defined(__GNUC__) || defined(__clang__)
                __builtin_prefetch(&nodes[n.children[0]], 0, 0);
                __builtin_prefetch(&nodes[n.children[1]], 0, 0);
#endif
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
    std::uint32_t* out_min_observed,
    double* packed_d0,
    const double* D0_data)
{
    if (N == 0 || Ndim == 0) {
        if (out_max_observed) *out_max_observed = 0;
        if (out_min_observed) *out_min_observed = 0;
        (void)walls;
        return;
    }
    // kdtree Cycle 4: probe bracket — master-thread wall covers all phases.
    const bool kd_probe_on = rcp_cycle_probe_enabled();
    const std::uint64_t kd_tsc_fn_start = kd_probe_on ? rcp_read_cycles() : 0;
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
        auto& cached_max_idx = g_dstat_max_idx;
        auto& cached_min_idx = g_dstat_min_idx;
        auto& cached_med_idx = g_dstat_med_idx;
        auto& cached_N = g_dstat_N;
        auto& cached_D0 = g_dstat_D0; auto& cached_Dmid = g_dstat_Dmid;
        auto& cached_valid = g_dstat_valid;
        // Cycle 21d: deterministic invalidation on permutation. The ratio
        // check below stays as belt-and-braces, but it can false-negative
        // under size ties (see g_rcp_perm_epoch comment).
        auto& cached_epoch = g_dstat_epoch;
        bool hit = false;
        if (cached_valid && cached_N == N && cached_D0 != 0.0 &&
            cached_Dmid != 0.0 && cached_epoch == g_rcp_perm_epoch) {
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
            cached_epoch = g_rcp_perm_epoch;
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
    t_shell *= g_rcp_shell_scale_dyn;   // Cycle 21k adaptive shell

    // Precompute inv_box for the MIC math — same as get_pairs_nd_3 so that
    // borderline pairs (d ≈ r_c) classify identically.
    double inv_box_kd[KD_MAX_NDIM] = {0.0};
    for (std::size_t d = 0; d < Ndim; ++d) inv_box_kd[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;  // [235] minimum image only along periodic dimensions: 0 along a wall -> no wrap

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
        g_kd_setup_cycles.fetch_add(rcp_read_cycles() - kd_tsc_setup_start,
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
    // Cycle 21k: while the adaptive shell is widened, force the full-
    // rebuild path. The patch path's fixed Dmin/2 query-radius inflation
    // assumes theta <= Dmin/2; widened budgets (theta * s_dyn) can exceed
    // that, and the widened phase is the high-refresh early phase where
    // full rebuild wins anyway.
    const bool need_full_rebuild =
        g_rcp_ckpt_deterministic ||          // [log 044/B4] no tree state across steps
        (g_kd_tree_N != N) ||
        (g_kd_patches_since_rebuild >= KD_MAX_PATCHES_BETWEEN_REBUILDS) ||
        (refresh_count_total * KD_PATCH_REFRESH_THRESHOLD_FRAC > N) ||
        (g_rcp_shell_scale_dyn > 1.0);

    const std::uint64_t kd_tsc_build_start = kd_probe_on ? rcp_read_cycles() : 0;
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
        kd_probe_on ? (rcp_read_cycles() - kd_tsc_build_start) : 0;
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
    //
    // Cycle 20 L2: the bbox-prune safety bound requires inflation >= theta
    // (an unrefreshed particle's true position can lag its leaf bbox by up
    // to the refresh trigger theta). theta = RCP_REFRESH_TRIGGER_SCALE *
    // Dmin/4; the historical Dmin/2 covers scale <= 2 exactly, so keep
    // Dmin/2 as the floor and add 25% margin above theta when the scale
    // pushes past 1.6.
    static const double kd_trigger_scale = []() {
        const char* s = std::getenv("RCP_REFRESH_TRIGGER_SCALE");
        if (!s || !s[0]) return 1.0;
        const double v = std::atof(s);
        return (v > 0.0 && v <= 2.0) ? v : 1.0;
    }();
    const double kd_radius_inflation =
        need_full_rebuild ? 0.0
                          : std::max(Dmin * 0.5,
                                     Dmin * 0.25 * kd_trigger_scale * 1.25);

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
    std::string overflow_msg;   // set under the overflow criticals below

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
    const std::uint64_t kd_tsc_qregion_start = kd_probe_on ? rcp_read_cycles() : 0;
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
        const std::uint64_t tls_kd_tsc_thread_start = kd_probe_on ? rcp_read_cycles() : 0;

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
            // Cycle 21i: the image radius must cover the FULL query reach
            // including the patch-path bbox inflation — the queries accept
            // through stale bboxes up to t_shell + inflation, so a
            // boundary-straddling pair reachable only via the inflation
            // margin needs its image enumerated too. The missing inflation
            // here was caught by the per-step missed-pair sentinel as rare
            // boundary-pair misses in 2D (interaction range is a larger
            // fraction of the box there). Strictly conservative fix.
            const std::uint64_t tsc_img0 = kd_probe_on ? rcp_read_cycles() : 0;
            const double* xi = x + i * Ndim;
            const ImageSet img = compute_image_set(
                xi, r_c_max + kd_radius_inflation, box.data(), Ndim);
            if (kd_probe_on) tls_kd_image += rcp_read_cycles() - tsc_img0;

            // Iterate the cartesian product of per-axis image offsets.
            // We unroll up to Ndim=8 here in a generic way using indices.
            std::size_t idx[KD_MAX_NDIM] = {0};
            for (std::uint32_t img_k = 0; img_k < img.total_images; ++img_k) {
                // Build xq for this image
                const std::uint64_t tsc_img1 = kd_probe_on ? rcp_read_cycles() : 0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    const std::int8_t off = img.offsets[d][idx[d]];
                    xq[d] = xi[d] + static_cast<double>(off) * box[d];
                }
                if (kd_probe_on) tls_kd_image += rcp_read_cycles() - tsc_img1;

                candidates.clear();
                const std::uint64_t tsc_q0 = kd_probe_on ? rcp_read_cycles() : 0;
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
                    tls_kd_query += rcp_read_cycles() - tsc_q0;
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
                    const std::uint64_t tsc_m0 = kd_probe_on ? rcp_read_cycles() : 0;
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
                    if (kd_probe_on) tls_kd_mic += rcp_read_cycles() - tsc_m0;
                    if (reject) continue;
                    if (kd_probe_on) ++tls_kd_acc;

                    const std::uint64_t tsc_w0 = kd_probe_on ? rcp_read_cycles() : 0;
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
                                warningIssued = true;
                                overflow_msg = rcp_nb_overflow_msg(
                                    "kdtree own-write", i, cap_i);
                                if (std::getenv("RCP_DEBUG_LAYOUT")) {
                                    std::cerr << "[ovf] i=" << i
                                              << " cap=" << cap_i
                                              << " Di=" << Di
                                              << " Dmax=" << Dmax
                                              << " t_shell=" << t_shell
                                              << " infl=" << kd_radius_inflation
                                              << " full_rebuild=" << need_full_rebuild
                                              << " kd_R=" << kd_R << "\n  slots:";
                                    const std::uint32_t* sl = pairs_data + pair_offsets[i];
                                    for (std::uint32_t kk = 1; kk <= std::min<std::uint32_t>(sl[0], 24); ++kk)
                                        std::cerr << " " << sl[kk];
                                    std::cerr << " ...\n";
                                }
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
                    if (kd_probe_on) tls_kd_write += rcp_read_cycles() - tsc_w0;
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
            tls_kd_total = rcp_read_cycles() - tls_kd_tsc_thread_start;
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
            // Dedup BEFORE the capacity check, matching the three sibling
            // append sites (kdtree own-write above, and both sort+walk
            // own/cross-write). Multi-image KD queries can enqueue the same
            // (j, i_writer) pair from several shifted images; the later
            // sort+unique only dedups AFTER this merge, so without a dedup
            // here transient duplicates consume slots and can force a FALSE
            // overflow in tight-capacity rows (e.g. small dense periodic
            // cells where cap clamps to N). Skipping known duplicates changes
            // nothing for normal runs (the row ends deduped either way) — it
            // only removes the spurious overflow. xwrite buffers stay tiny, so
            // the added scan cost is negligible.
            bool dupj = false;
            for (std::uint32_t kk = 1; kk <= cntj; ++kk) {
                if (slot_j[kk] == i_writer) { dupj = true; break; }
            }
            if (dupj) continue;
            if (cntj < cap_j) {
                slot_j[0] = cntj + 1;
                slot_j[cntj + 1] = i_writer;
                if (!dirty_kd[j]) {
                    kd_dirty_list.push_back(j);
                    dirty_kd[j] = 1;
                }
            } else if (!warningIssued) {
                warningIssued = true;
                overflow_msg = rcp_nb_overflow_msg(
                    "kdtree cross-write", j, cap_j);
            }
        }
    }

    // Cycle 23 probe (gated by RCP_DEBUG_LAYOUT): confirms the cross-write
    // buffers stay tiny (~1 MB) and reports pairs_data reserved size. Verified
    // xwrite is NOT the memory hog; the resident growth was pairs_data faulted
    // across shifting reorder layouts, fixed via LazyPairs::release_resident().
    static const bool kd_xwrite_probe = []() {
        const char* s = std::getenv("RCP_DEBUG_LAYOUT");
        return s && s[0] && s[0] != '0';
    }();
    if (kd_xwrite_probe) {
        static std::uint64_t s_peak_live = 0;
        std::uint64_t live = 0, cap = 0;
        for (int t = 0; t < kd_max_threads; ++t) {
            live += g_kd_xwrite_buffers[t].size();
            cap  += g_kd_xwrite_buffers[t].capacity();
        }
        if (live > s_peak_live) {
            s_peak_live = live;
            std::fprintf(stderr,
                "[xwrite] kd_R=%zu live=%.3fGB cap=%.3fGB pairs_data_cap=%.3fGB\n",
                kd_R, live * 8.0 / 1e9, cap * 8.0 / 1e9,
                static_cast<double>(pair_offsets[N]) * 4.0 / 1e9);
        }
    }

    // Cycle 21N: surface any recorded overflow (parallel own-write or the
    // serial cross-write merge above) as an exception now that all OpenMP
    // regions and locks are released.
    if (warningIssued) throw std::runtime_error(overflow_msg);

    // Cycle 7 measurement: record query-region wall per regime.
    if (kd_probe_on) {
        const std::uint64_t qregion = rcp_read_cycles() - kd_tsc_qregion_start;
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
    const std::uint64_t kd_tsc_sort_start = kd_probe_on ? rcp_read_cycles() : 0;
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
            cnt = slot[0];
        }
        // Cycle 21k: refill the packed-D0 mirror. The dirty list covers
        // every row this call mutated (refresh rebuilds at the accept-write
        // sites AND cross-write appends, which push their target onto the
        // list), so the mirror can never go stale.
        if (packed_d0) {
            double* pk = packed_d0 + pair_offsets[i];
            for (std::uint32_t k = 1; k <= cnt; ++k)
                pk[k] = D0_data[slot[k]];
        }
    }
    if (kd_probe_on) {
        g_kd_sort_cycles.fetch_add(rcp_read_cycles() - kd_tsc_sort_start,
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
        g_kd_total_cycles.fetch_add(rcp_read_cycles() - kd_tsc_fn_start,
                                    std::memory_order_relaxed);
        g_kd_calls.fetch_add(1, std::memory_order_relaxed);
    }
}

// ============================================================================
// [227] N-dimensional k-d tree pair search, used only for Ndim >= 4.
//
// 2D and 3D keep their own optimised backends untouched. For Ndim >= 4 the only
// path used to be sort+walk (sort along x, scan an x-window), whose candidate count
// per particle grows like N * (reach / L): in 4-6D that is a large fraction of N.
// This backend replaces only the CANDIDATE GENERATION: a dimension-generic k-d tree
// (double bounding boxes of Ndim entries, per-node max diameter for tight pruning,
// median splits cycling through the axes). Every candidate is then judged by exactly
// the sort+walk acceptance test (same shell t, r_c = Di/2 + t + float(Dj/2), same
// dx-first squared distance, minimum image on every axis) and stored by exactly the
// sort+walk rules (pair kept in the lower-index particle's row; refreshed rows rebuilt;
// cross-writes into unrefreshed rows; every row sorted), so the neighbour lists are
// identical to sort+walk's and runs are bit-for-bit the same, only faster. The tree is
// rebuilt on every call (O(N log N), parallel); the queries run in parallel over the
// refreshed particles with per-thread deferred cross-writes, so results do not depend
// on the thread count. RCP_PAIRS_BACKEND=sortwalk selects the old path; =shadow runs
// both and compares the pair sets on every call.
// ============================================================================
namespace {
struct KdTreeND {
    static constexpr std::uint32_t LEAF = 32;
    static constexpr std::uint32_t PAR_MIN = 4096;   // below this subtree size, build serially
    std::size_t nd = 0;
    std::vector<std::uint32_t> perm, child0, child1, first, count;
    std::vector<double> lo, hi, maxD;                // lo/hi: nd entries per node
    std::atomic<std::uint32_t> next{0};

    void build(const double* x, const double* D, std::size_t N, std::size_t Ndim) {
        nd = Ndim;
        perm.resize(N);
        for (std::size_t k = 0; k < N; ++k) perm[k] = static_cast<std::uint32_t>(k);
        const std::size_t upper = 4 * N / LEAF + 64;
        if (child0.size() < upper) {
            child0.resize(upper); child1.resize(upper); first.resize(upper); count.resize(upper);
            maxD.resize(upper);
        }
        if (lo.size() < upper * nd) { lo.resize(upper * nd); hi.resize(upper * nd); }
        next.store(1, std::memory_order_relaxed);
        built_n = N;
        if (N == 0) return;
        #pragma omp parallel
        {
            #pragma omp single
            build_rec(x, D, 0, 0, static_cast<std::uint32_t>(N), 0);
        }
        node_count = next.load(std::memory_order_relaxed);
        leaves.clear();
        for (std::uint32_t n = 0; n < node_count; ++n)
            if (child0[n] == KD_NULL) leaves.push_back(n);
    }

    // Keep the topology, recompute every box and max diameter from the current positions and
    // diameters: O(N Ndim). Boxes stay exact, so queries remain a superset of the exact test;
    // only pruning efficiency degrades as particles drift, which a periodic rebuild restores.
    // Children always have larger indices than their parent, so a reverse sweep is bottom-up.
    void refit(const double* x, const double* D) {
        #pragma omp parallel for schedule(static)
        for (std::size_t k = 0; k < leaves.size(); ++k) {
            const std::uint32_t node = leaves[k];
            double* nlo = &lo[std::size_t(node) * nd];
            double* nhi = &hi[std::size_t(node) * nd];
            for (std::size_t d = 0; d < nd; ++d) {
                nlo[d] = std::numeric_limits<double>::infinity();
                nhi[d] = -std::numeric_limits<double>::infinity();
            }
            double mx = 0.0;
            for (std::uint32_t q = first[node]; q < first[node] + count[node]; ++q) {
                const std::uint32_t p = perm[q];
                const double* xp = x + std::size_t(p) * nd;
                for (std::size_t d = 0; d < nd; ++d) {
                    if (xp[d] < nlo[d]) nlo[d] = xp[d];
                    if (xp[d] > nhi[d]) nhi[d] = xp[d];
                }
                if (D[p] > mx) mx = D[p];
            }
            maxD[node] = mx;
        }
        for (std::uint32_t n = node_count; n-- > 0;) {
            if (child0[n] == KD_NULL) continue;
            const std::uint32_t l = child0[n], r = child1[n];
            for (std::size_t d = 0; d < nd; ++d) {
                lo[std::size_t(n) * nd + d] = std::min(lo[std::size_t(l) * nd + d], lo[std::size_t(r) * nd + d]);
                hi[std::size_t(n) * nd + d] = std::max(hi[std::size_t(l) * nd + d], hi[std::size_t(r) * nd + d]);
            }
            maxD[n] = std::max(maxD[l], maxD[r]);
        }
    }

    std::size_t built_n = 0;
    std::uint32_t node_count = 0;
    std::vector<std::uint32_t> leaves;

    // Candidates for particle position xq: every particle p whose minimum-image distance is
    // <= (Di + D[p]) / 2 + reach. One traversal covers every periodic image: the distance from xq
    // to a node's box is the minimum over the box and its copies shifted by +-L on each axis (a
    // lower bound on the minimum-image distance to any particle inside). A superset of the exact
    // test applied by the caller.
    void query(const double* xq, double Di, double reach, const double* x, const double* D,
               const double* L, const double* invL, std::vector<std::uint32_t>& out) const {
        std::uint32_t stack[128];
        int top = 0;
        stack[top++] = 0;
        while (top > 0) {
            const std::uint32_t n = stack[--top];
            const double* nlo = &lo[std::size_t(n) * nd];
            const double* nhi = &hi[std::size_t(n) * nd];
            double dbb2 = 0.0;
            for (std::size_t d = 0; d < nd; ++d) {
                const double direct = std::max(0.0, std::max(nlo[d] - xq[d], xq[d] - nhi[d]));
                const double wrapped = std::min(xq[d] + L[d] - nhi[d], nlo[d] - xq[d] + L[d]);
                const double diff = std::max(0.0, std::min(direct, wrapped));
                dbb2 += diff * diff;
            }
            const double rn = (Di + maxD[n]) * 0.5 + reach;
            if (dbb2 > rn * rn) continue;
            if (child0[n] == KD_NULL) {
                const std::uint32_t f = first[n], c = count[n];
                for (std::uint32_t k = 0; k < c; ++k) {
                    const std::uint32_t p = perm[f + k];
                    const double* xp = x + std::size_t(p) * nd;
                    double d2 = 0.0;
                    for (std::size_t d = 0; d < nd; ++d) {
                        double dq = xp[d] - xq[d];
                        dq -= std::nearbyint(dq * invL[d]) * L[d];
                        d2 += dq * dq;
                    }
                    const double rp = (Di + D[p]) * 0.5 + reach;
                    if (d2 <= rp * rp) out.push_back(p);
                }
                continue;
            }
            stack[top++] = child0[n];
            stack[top++] = child1[n];
        }
    }

private:
    void build_rec(const double* x, const double* D, std::uint32_t node,
                   std::uint32_t b, std::uint32_t e, std::uint32_t depth) {
        const std::uint32_t cnt = e - b;
        double* nlo = &lo[std::size_t(node) * nd];
        double* nhi = &hi[std::size_t(node) * nd];
        if (cnt <= LEAF) {
            for (std::size_t d = 0; d < nd; ++d) {
                nlo[d] = std::numeric_limits<double>::infinity();
                nhi[d] = -std::numeric_limits<double>::infinity();
            }
            double mx = 0.0;
            for (std::uint32_t k = b; k < e; ++k) {
                const std::uint32_t p = perm[k];
                const double* xp = x + std::size_t(p) * nd;
                for (std::size_t d = 0; d < nd; ++d) {
                    if (xp[d] < nlo[d]) nlo[d] = xp[d];
                    if (xp[d] > nhi[d]) nhi[d] = xp[d];
                }
                if (D[p] > mx) mx = D[p];
            }
            maxD[node] = mx; first[node] = b; count[node] = cnt;
            child0[node] = KD_NULL; child1[node] = KD_NULL;
            return;
        }
        const std::size_t sd = depth % nd;
        const std::uint32_t mid = b + cnt / 2;
        const std::size_t ndim = nd;
        std::nth_element(perm.begin() + b, perm.begin() + mid, perm.begin() + e,
                         [x, ndim, sd](std::uint32_t a, std::uint32_t c) {
                             return x[std::size_t(a) * ndim + sd] < x[std::size_t(c) * ndim + sd];
                         });
        const std::uint32_t l = next.fetch_add(2, std::memory_order_relaxed);
        const std::uint32_t r = l + 1;
        child0[node] = l; child1[node] = r; count[node] = 0;
        if (cnt >= PAR_MIN) {
            #pragma omp task firstprivate(l, b, mid, depth) shared(x, D)
            build_rec(x, D, l, b, mid, depth + 1);
            #pragma omp task firstprivate(r, mid, e, depth) shared(x, D)
            build_rec(x, D, r, mid, e, depth + 1);
            #pragma omp taskwait
        } else {
            build_rec(x, D, l, b, mid, depth + 1);
            build_rec(x, D, r, mid, e, depth + 1);
        }
        const double* llo = &lo[std::size_t(l) * nd]; const double* lhi = &hi[std::size_t(l) * nd];
        const double* rlo = &lo[std::size_t(r) * nd]; const double* rhi = &hi[std::size_t(r) * nd];
        for (std::size_t d = 0; d < nd; ++d) {
            nlo[d] = std::min(llo[d], rlo[d]);
            nhi[d] = std::max(lhi[d], rhi[d]);
        }
        maxD[node] = std::max(maxD[l], maxD[r]);
    }
};

// RCP_KDND_STATS=1: counters for the cost model of this backend, printed to stderr at process exit.
// Off by default; when off nothing is counted or timed.
struct KdNdStats {
    bool on = false;
    std::uint64_t calls = 0, refreshed = 0, candidates = 0, accepted = 0, n = 0, ndim = 0, rebuilds = 0;
    double build_s = 0.0, query_s = 0.0;
    KdNdStats() {
        const char* s = std::getenv("RCP_KDND_STATS");
        on = s && s[0] && s[0] != '0';
    }
    ~KdNdStats() {
        if (!on || calls == 0) return;
        const double q = refreshed ? static_cast<double>(refreshed) : 1.0;
        std::fprintf(stderr,
                     "[kdnd-stats] Ndim=%llu N=%llu calls=%llu rebuilds=%llu refreshed_per_call=%.2f "
                     "candidates_per_query=%.2f accepted_per_query=%.2f build_ms_per_call=%.4f "
                     "query_us_per_particle=%.4f\n",
                     static_cast<unsigned long long>(ndim), static_cast<unsigned long long>(n),
                     static_cast<unsigned long long>(calls), static_cast<unsigned long long>(rebuilds),
                     refreshed / static_cast<double>(calls),
                     candidates / q, accepted / q, 1e3 * build_s / calls, 1e6 * query_s / q);
    }
};
static KdNdStats g_kdnd_stats;
}  // namespace

void get_pairs_kdtree_nd(
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
    std::uint32_t* out_min_observed,
    double* packed_d0,
    const double* D0_data)
{
    (void)walls;   // like sort+walk: minimum image on every axis (walls keep particles apart)
    if (N == 0) {
        if (out_max_observed) *out_max_observed = 0;
        if (out_min_observed) *out_min_observed = 0;
        ++g_pairs_version;
        return;
    }
    // Shell thickness exactly as sort+walk computes it.
    const double Dmax = *std::max_element(D.begin(), D.end());
    const double Dmin = *std::min_element(D.begin(), D.end());
    std::vector<double> Dsorted = D;
    std::nth_element(Dsorted.begin(), Dsorted.begin() + N / 2, Dsorted.end());
    const double median = Dsorted[N / 2];
    double t = std::min(std::max(median * 1.02, Dmin * 2.25), Dmax) * 0.5;
    t *= g_rcp_shell_scale_dyn;   // Cycle 21k adaptive shell
    std::vector<double> inv_box(Ndim);
    for (std::size_t d = 0; d < Ndim; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;  // [235] minimum image only along periodic dimensions: 0 along a wall -> no wrap
    // The tree proposes a superset: its Euclidean, double-precision test gets a small margin
    // over the exact test (float rounding of Dj/2, image-shift rounding).
    const double margin = 1e-6 * Dmax;

    std::vector<std::uint32_t> refreshed;
    refreshed.reserve(N);
    for (std::size_t k = 0; k < N; ++k)
        if (refresh[k] == 1) refreshed.push_back(static_cast<std::uint32_t>(k));
    const std::size_t R = refreshed.size();

    // Tree maintenance: refit the existing topology (O(N Ndim)) on most calls; rebuild (O(N log N))
    // when the particle set or dimension changed, every 32 calls, when more than 1/20 of the
    // particles were refreshed (early, mobile regime), or while the adaptive shell is widened.
    // The pair lists do not depend on which is done: the exact test decides membership.
    const bool stats = g_kdnd_stats.on;
    const auto tb0 = std::chrono::steady_clock::now();
    static KdTreeND tree;
    static std::size_t calls_since_rebuild = 0;
    static std::uint64_t built_epoch = ~0ULL;   // particle renumbering (spatial reorder) scatters leaves
    const bool rebuild = tree.built_n != N || tree.nd != Ndim || calls_since_rebuild >= 32 ||
                         R * 20 > N || g_rcp_shell_scale_dyn > 1.0 || built_epoch != g_rcp_perm_epoch;
    if (rebuild) {
        tree.build(x, D.data(), N, Ndim);
        calls_since_rebuild = 0;
        built_epoch = g_rcp_perm_epoch;
    }
    else { tree.refit(x, D.data()); ++calls_since_rebuild; }
    const auto tb1 = std::chrono::steady_clock::now();
    if (stats) { if (rebuild) ++g_kdnd_stats.rebuilds; }

    const int nthreads = omp_get_max_threads();
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> xwrite(nthreads);
    bool overflow = false;
    std::string overflow_msg;

    std::uint64_t cand_total = 0, acc_total = 0;
    #pragma omp parallel reduction(+:cand_total, acc_total)
    {
        auto& xbuf = xwrite[omp_get_thread_num()];
        std::vector<std::uint32_t> cand;
        cand.reserve(512);
        #pragma omp for schedule(dynamic, 32)
        for (std::size_t r = 0; r < R; ++r) {
            const std::size_t i = refreshed[r];
            std::uint32_t* slot_i = pairs_data + pair_offsets[i];
            slot_i[0] = 0;
            const std::uint32_t cap_i = max_neighbors_per_particle[i];
            const double Di = D[i];
            const double Di_half_plus_t = Di * 0.5 + t;
            const double* xi = x + i * Ndim;
            cand.clear();
            tree.query(xi, Di, t + margin, x, D.data(), box.data(), inv_box.data(), cand);
            cand_total += cand.size();
            {
                for (const std::uint32_t jj : cand) {
                    const std::size_t j = jj;
                    if (j == i) continue;
                    // the exact sort+walk acceptance test
                    double dx = x[j * Ndim] - xi[0];
                    dx -= std::nearbyint(dx * inv_box[0]) * box[0];
                    const double r_c = Di_half_plus_t + static_cast<double>(static_cast<float>(D[j] * 0.5));
                    if (!(std::fabs(dx) < r_c)) continue;
                    double d2 = dx * dx;
                    for (std::size_t d = 1; d < Ndim; ++d) {
                        double dz = x[j * Ndim + d] - xi[d];
                        dz -= std::nearbyint(dz * inv_box[d]) * box[d];
                        d2 += dz * dz;
                    }
                    if (d2 >= r_c * r_c) continue;
                    ++acc_total;
                    if (j > i) {
                        std::uint32_t& cnt = slot_i[0];
                        bool dup = false;   // defensive: a row never holds the same j twice
                        for (std::uint32_t kk = 1; kk <= cnt; ++kk)
                            if (slot_i[kk] == jj) { dup = true; break; }
                        if (dup) continue;
                        if (cnt < cap_i) {
                            ++cnt;
                            slot_i[cnt] = jj;
                        } else {
                            #pragma omp critical(warn_neighbor_overflow_kdnd)
                            if (!overflow) {
                                overflow = true;
                                overflow_msg = rcp_nb_overflow_msg("kdtree-nd own-write", i, cap_i);
                            }
                        }
                    } else if (refresh[j] == 0) {
                        xbuf.emplace_back(jj, static_cast<std::uint32_t>(i));
                    }
                    // else: j < i and refreshed -> j's own query stores the pair
                }
            }
        }
    }
    if (stats) {
        const auto tq1 = std::chrono::steady_clock::now();
        g_kdnd_stats.calls += 1; g_kdnd_stats.refreshed += R; g_kdnd_stats.n = N; g_kdnd_stats.ndim = Ndim;
        g_kdnd_stats.candidates += cand_total; g_kdnd_stats.accepted += acc_total;
        g_kdnd_stats.build_s += std::chrono::duration<double>(tb1 - tb0).count();
        g_kdnd_stats.query_s += std::chrono::duration<double>(tq1 - tb1).count();
    }
    // rows touched this call: the refreshed ones (rebuilt) and the cross-write recipients. Only
    // these can be unsorted; every other row is unchanged since its last sort (as in the 3D kd-tree).
    std::vector<std::uint32_t> dirty(refreshed);
    static std::vector<std::uint8_t> is_dirty;
    if (is_dirty.size() != N) is_dirty.assign(N, 0);
    for (const std::uint32_t k : refreshed) is_dirty[k] = 1;
    // cross-writes into unrefreshed rows, serially, with the same dedup and capacity rule
    for (int th = 0; th < nthreads && !overflow; ++th) {
        for (const auto& jp : xwrite[th]) {
            const std::uint32_t j = jp.first, iw = jp.second;
            std::uint32_t* slot_j = pairs_data + pair_offsets[j];
            const std::uint32_t cntj = slot_j[0];
            bool dup = false;
            for (std::uint32_t kk = 1; kk <= cntj; ++kk)
                if (slot_j[kk] == iw) { dup = true; break; }
            if (dup) continue;
            if (cntj < max_neighbors_per_particle[j]) {
                slot_j[0] = cntj + 1;
                slot_j[cntj + 1] = iw;
                if (!is_dirty[j]) { is_dirty[j] = 1; dirty.push_back(j); }
            } else {
                overflow = true;
                overflow_msg = rcp_nb_overflow_msg("kdtree-nd cross-write", j, max_neighbors_per_particle[j]);
                break;
            }
        }
    }
    for (const std::uint32_t k : dirty) is_dirty[k] = 0;   // reset for the next call (also on overflow)
    if (overflow) throw std::runtime_error(overflow_msg);

    // canonical form, as sort+walk: every row sorted (only touched rows can be out of order)
    #pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t r = 0; r < dirty.size(); ++r) {
        const std::size_t i = dirty[r];
        std::uint32_t* slot = pairs_data + pair_offsets[i];
        const std::uint32_t cnt = slot[0];
        if (cnt > 1) std::sort(slot + 1, slot + 1 + cnt);
        if (packed_d0) {
            double* pk = packed_d0 + pair_offsets[i];
            for (std::uint32_t k = 1; k <= cnt; ++k) pk[k] = D0_data[slot[k]];
        }
    }
    if (out_max_observed || out_min_observed) {
        std::uint32_t lmax = 0, lmin = std::numeric_limits<std::uint32_t>::max();
        #pragma omp parallel for schedule(static) reduction(max:lmax) reduction(min:lmin)
        for (std::size_t i = 0; i < N; ++i) {
            const std::uint32_t c = pairs_data[pair_offsets[i]];
            if (c > lmax) lmax = c;
            if (c < lmin) lmin = c;
        }
        if (out_max_observed) *out_max_observed = lmax;
        if (out_min_observed) *out_min_observed = lmin;
    }
    ++g_pairs_version;
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
    std::uint32_t* out_min_observed,
    double* packed_d0,
    const double* D0_data)
{
    PairsBackend backend = rcp_pairs_backend();
    // Cycle 21h: the 2D/3D kdtree is structurally Ndim <= 3 (H18 float bbox[3]).
    // [227] Ndim >= 4 has its own dimension-generic k-d tree (get_pairs_kdtree_nd), which
    // reproduces sort+walk's neighbour lists exactly; RCP_PAIRS_BACKEND=sortwalk selects the
    // old path and =shadow runs both and compares every call.
    if (Ndim > 3) {
        if (backend == PairsBackend::SortWalk) {
            get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                           max_neighbors_per_particle, pair_offsets, pairs_data,
                           out_max_observed, out_min_observed, packed_d0, D0_data);
            return;
        }
        if (backend == PairsBackend::KdTree) {
            get_pairs_kdtree_nd(N, Ndim, x, D, box, walls, refresh,
                                max_neighbors_per_particle, pair_offsets, pairs_data,
                                out_max_observed, out_min_observed, packed_d0, D0_data);
            return;
        }
        const std::uint64_t total_nd = pair_offsets[N];
        static thread_local std::vector<std::uint32_t> shadow_nd;
        if (shadow_nd.size() != total_nd) shadow_nd.resize(total_nd);
        std::memcpy(shadow_nd.data(), pairs_data, total_nd * sizeof(std::uint32_t));
        std::uint32_t smax = 0, smin = std::numeric_limits<std::uint32_t>::max();
        get_pairs_kdtree_nd(N, Ndim, x, D, box, walls, refresh,
                            max_neighbors_per_particle, pair_offsets, shadow_nd.data(),
                            &smax, &smin, nullptr, nullptr);
        get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                       max_neighbors_per_particle, pair_offsets, pairs_data,
                       out_max_observed, out_min_observed, packed_d0, D0_data);
        compare_pairs_sets(N, pair_offsets, pairs_data, shadow_nd.data());
        return;
    }
    if (backend == PairsBackend::SortWalk) {
        get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                       max_neighbors_per_particle, pair_offsets, pairs_data,
                       out_max_observed, out_min_observed,
                       packed_d0, D0_data);
        return;
    }
    if (backend == PairsBackend::KdTree) {
        get_pairs_kdtree(N, Ndim, x, D, box, walls, refresh,
                         max_neighbors_per_particle, pair_offsets, pairs_data,
                         out_max_observed, out_min_observed,
                         packed_d0, D0_data);
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
                     &shadow_max, &shadow_min,
                     nullptr, nullptr);   // shadow run never touches packed
    get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh,
                   max_neighbors_per_particle, pair_offsets, pairs_data,
                   out_max_observed, out_min_observed,
                   packed_d0, D0_data);
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
    std::vector<std::size_t>& z,
    const double* packed_d0,
    double kappa,
    // Cascade (additive): per-particle "kappa radius" = fixed?0:D[i]/2, and the
    // active-only diameter sum. nullptr / negative -> standard path (byte-identical):
    // the masked dkappa terms below reduce EXACTLY to the originals.
    const double* rk_data,
    double sumD_active,
    SwapEdgeSink* swap_sink)
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
    // Cycle 21e: per-thread touched-particle ranges from the PREVIOUS call.
    // Each thread writes its own static-schedule i-chunk plus cross-writes
    // to j's, which under Morton reorder are index-local — so the touched
    // range per thread is far smaller than N. Zero (here) and reduce
    // (below) only the touched range; the rest of the buffer is known-zero
    // from the previous call's zeroing (it was inside that call's range or
    // never written since the initial full zero). Falls back to full range
    // whenever the buffer is (re)allocated or the thread count changes.
    constexpr int RCP_MAX_THREADS = 64;
    static std::size_t g_fl_lo[RCP_MAX_THREADS];
    static std::size_t g_fl_hi[RCP_MAX_THREADS];   // exclusive
    static int g_fl_threads = -1;
    static std::size_t g_fl_N = 0;
    const bool ranges_valid =
        (g_fl_threads == nthreads) && (g_fl_N == N) &&
        nthreads <= RCP_MAX_THREADS;
    #pragma omp parallel num_threads(nthreads)
    {
        const int tid = omp_get_thread_num();
        double* Ft = F_base + tid * F_stride;
        std::size_t* zt = z_base + tid * z_stride;
        if (ranges_valid) {
            const std::size_t lo = g_fl_lo[tid];
            const std::size_t hi = g_fl_hi[tid];
            if (lo < hi) {
                std::fill(Ft + lo * Ndim, Ft + hi * Ndim, 0.0);
                std::fill(zt + lo, zt + hi, std::size_t{0});
            }
        } else {
            std::fill(Ft, Ft + F_stride, 0.0);
            std::fill(zt, zt + z_stride, std::size_t{0});
        }
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
    for (std::size_t d = 0; d < Ndim; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box_data[d] : 0.0;  // [235] minimum image only along periodic dimensions: 0 along a wall -> no wrap

    g_t_forces_loop.begin();
    const bool probe = rcp_cycle_probe_enabled();
    // Cycle 15c MIC-cost probe. Each extra unit injects one extra
    // `floor(delta*inv_box + 0.5) * box` computation per dim, scaled by 1e-16
    // before adding to delta — runtime nonzero so the compiler can't elide,
    // but small enough to be ~1 ulp noise on a delta ~ O(1). Slope of
    // cyc/cand vs extra_mic = per-MIC cost.
    const int extra_mic = rcp_extra_mic_count();
    const double mic_scale = 1e-16;
    // Cycle 21: SIMD reject-test eligibility, resolved once per call.
    const bool use_simd3 =
#if defined(__AVX512F__)
        rcp_force_simd() && Ndim == 3 && extra_mic == 0;
#else
        false;
#endif
    // Cycle 21f: 2D variant (same structure, two coordinate gathers).
    const bool use_simd2 =
#if defined(__AVX512F__)
        rcp_force_simd() && Ndim == 2 && extra_mic == 0;
#else
        false;
#endif
    // Cycle 21k: packed-D0 contiguous loads replace the D[j] gather in the
    // SIMD reject paths. Dj = round(packed*kappa) is the same single
    // rounding as the D-update pass's D[j] = round(D0[j]*kappa), so accept
    // decisions and forces are byte-identical to the gather path.
    const bool use_packed = packed_d0 != nullptr;
    // Cycle 21L FP32 stub: per-call float shadow of x. The real engine
    // would fold these stores into the fused position pass; the rebuild
    // here costs ~2% of floop and is included in the measured number
    // (i.e. the stub UNDERSTATES the engine's gain slightly).
    const bool use_fp32_stub = rcp_fp32_stub() && use_simd3;
    static std::vector<float> g_x_shadow;
    if (use_fp32_stub) {
        const std::size_t total = N * Ndim;
        if (g_x_shadow.size() != total) g_x_shadow.resize(total);
        float* xs = g_x_shadow.data();
        #pragma omp parallel for schedule(static)
        for (std::size_t k = 0; k < total; ++k)
            xs[k] = static_cast<float>(x_data[k]);
    }
    const float* x_shadow = use_fp32_stub ? g_x_shadow.data() : nullptr;
    // SWAP near-pair collection (only when a sink is supplied, i.e. on sweep steps).
    const double swap_margin = swap_sink ? swap_sink->margin : 0.0;
    if (swap_sink) {
        // capacity = every pair-list entry (worst case all near)
        std::size_t total_entries = 0;
        if (swap_sink->want_edges)
            for (std::size_t i = 0; i < N; ++i) total_entries += pairs_data[pair_offsets[i]];
        swap_sink->reset(nthreads, total_entries, N);
    }
    #pragma omp parallel default(none) \
        shared(x_data, pairs_data, pair_offsets, D_data, box_data, inv_box, \
               F_data, z_data, N, Ndim, K, F_stride, z_stride, \
               g_floop_cycles, g_floop_candidates, g_floop_overlaps, \
               g_floop_test_cycles, g_floop_overlap_cycles, probe, \
               extra_mic, mic_scale, use_simd3, use_simd2, \
               g_fl_lo, g_fl_hi, packed_d0, kappa, use_packed, \
               use_fp32_stub, x_shadow, rk_data, swap_sink, swap_margin) \
        reduction(+:U_red, Lc_red, dkappa_red, Fmean_red, count_red) \
        reduction(max:max_min_red)
    {
        const int tid = omp_get_thread_num();
        double*      F_t = F_data + tid * F_stride;
        std::size_t* z_t = z_data + tid * z_stride;
        SwapEdge* swap_edges = (swap_sink && swap_sink->want_edges) ? swap_sink->per_thread[tid].data() : nullptr;
        std::size_t swap_n = 0;
        double* swap_g = (swap_sink && swap_sink->want_gsum)
            ? swap_sink->gsum_thread.data() + static_cast<std::size_t>(tid) * swap_sink->gstride : nullptr;
        std::uint64_t tls_cycles = 0;
        std::uint64_t tls_candidates = 0;
        std::uint64_t tls_overlaps = 0;
        // Cycle 21e: track this thread's touched-particle index range
        // (i writes + j cross-writes) for range-limited zero/reduce.
        std::size_t tls_lo = N, tls_hi = 0;
        // Cycle 15c: split per-candidate cycles into "test path" (always
        // executed: load xj, MIC wrap, d2, compare) and "overlap path"
        // (sqrt + F/z buffer writes). Lets us see which side dominates
        // and whether the 38 cyc/cand surplus is in test or overlap.
        std::uint64_t tls_cyc_test = 0;
        std::uint64_t tls_cyc_overlap = 0;
        std::uint64_t tsc_start = probe ? rcp_read_cycles() : 0;

        #pragma omp for schedule(static)
        for (std::size_t i = 0; i < N; ++i) {
            const std::uint32_t* pi = pairs_data + pair_offsets[i];
            const double* xi = x_data + i * Ndim;
            const double Di = D_data[i];
            std::uint32_t numNbr = pi[0];
            // Cycle 6: pairs[i] now stores only neighbors j > i (asymmetric
            // list). The `if (j <= i) continue;` filter is no longer needed —
            // every entry contributes work.
            std::uint32_t jdx = 1;
#if defined(__AVX512F__)
            if (use_simd3) {
                // Cycle 21: 8-wide AVX-512 batching of the reject test.
                // Lane math mirrors the scalar path op-for-op (fnmadd for
                // the MIC subtract, fmadd for the d2 accumulation — the
                // contractions gcc emits for the scalar source) and lanes
                // are consumed in jdx order, so accept decisions and
                // accumulation order match the scalar path. Accepted pairs
                // run the identical scalar physics body.
                const __m512d xi0 = _mm512_set1_pd(xi[0]);
                const __m512d xi1 = _mm512_set1_pd(xi[1]);
                const __m512d xi2 = _mm512_set1_pd(xi[2]);
                const __m512d Di_v = _mm512_set1_pd(Di);
                const __m512d half_v = _mm512_set1_pd(0.5);
                const __m512d kap_v = _mm512_set1_pd(kappa);
                const double* pk_i =
                    use_packed ? packed_d0 + pair_offsets[i] : nullptr;
                const __m512d ib0 = _mm512_set1_pd(inv_box[0]);
                const __m512d ib1 = _mm512_set1_pd(inv_box[1]);
                const __m512d ib2 = _mm512_set1_pd(inv_box[2]);
                const __m512d bx0 = _mm512_set1_pd(box_data[0]);
                const __m512d bx1 = _mm512_set1_pd(box_data[1]);
                const __m512d bx2 = _mm512_set1_pd(box_data[2]);
                const __m256i one3 = _mm256_set1_epi32(3);
                // Cycle 21b: masked final block instead of a scalar tail.
                // Lists are short on average (~10-30 entries), so a scalar
                // tail leaves a large candidate fraction unvectorized; one
                // masked iteration covers it. Masked-off lanes gather from
                // index 0 (safe: slot 0 of x/D always exists) and are
                // excluded from the accept mask.
                for (std::uint32_t blk = 1; blk <= numNbr; blk += 8) {
                    jdx = blk;
                    const std::uint32_t remaining = numNbr - blk + 1;
                    const __mmask8 lanemask =
                        remaining >= 8
                            ? static_cast<__mmask8>(0xFF)
                            : static_cast<__mmask8>((1u << remaining) - 1u);
                    const __m256i jv = _mm256_maskz_loadu_epi32(
                        lanemask, pi + blk);
                    const __m256i j3 = _mm256_mullo_epi32(jv, one3);
                    const __m512d zero_v = _mm512_setzero_pd();
                    __m512d xj0, xj1, xj2;
                    if (use_fp32_stub) {
                        // Cycle 21L stub: 8-element FLOAT gathers (half the
                        // bytes), converted to double for the unchanged
                        // test math. Measurement-only path.
                        const __m256 zf = _mm256_setzero_ps();
                        xj0 = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(
                            zf, lanemask, j3, x_shadow, 4));
                        xj1 = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(
                            zf, lanemask,
                            _mm256_add_epi32(j3, _mm256_set1_epi32(1)),
                            x_shadow, 4));
                        xj2 = _mm512_cvtps_pd(_mm256_mmask_i32gather_ps(
                            zf, lanemask,
                            _mm256_add_epi32(j3, _mm256_set1_epi32(2)),
                            x_shadow, 4));
                    } else {
                        xj0 = _mm512_mask_i32gather_pd(
                            zero_v, lanemask, j3, x_data, 8);
                        xj1 = _mm512_mask_i32gather_pd(
                            zero_v, lanemask,
                            _mm256_add_epi32(j3, _mm256_set1_epi32(1)), x_data, 8);
                        xj2 = _mm512_mask_i32gather_pd(
                            zero_v, lanemask,
                            _mm256_add_epi32(j3, _mm256_set1_epi32(2)), x_data, 8);
                    }
                    // Cycle 21k: contiguous packed-D0 load (masked lanes
                    // are not accessed, so the row tail is fault-safe)
                    // instead of the D[j] gather. Masked-off lanes read 0
                    // here vs D[0] in the gather — both are excluded from
                    // the accept mask by lanemask, so the difference is
                    // unobservable.
                    const __m512d Dj = use_packed
                        ? _mm512_mul_pd(
                              _mm512_maskz_loadu_pd(lanemask, pk_i + blk),
                              kap_v)
                        : _mm512_mask_i32gather_pd(
                              zero_v, lanemask, jv, D_data, 8);

                    __m512d d0 = _mm512_sub_pd(xj0, xi0);
                    __m512d t = _mm512_roundscale_pd(
                        _mm512_mul_pd(d0, ib0),
                        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
                    d0 = _mm512_fnmadd_pd(t, bx0, d0);

                    __m512d d1 = _mm512_sub_pd(xj1, xi1);
                    t = _mm512_roundscale_pd(
                        _mm512_mul_pd(d1, ib1),
                        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
                    d1 = _mm512_fnmadd_pd(t, bx1, d1);

                    __m512d d2dim = _mm512_sub_pd(xj2, xi2);
                    t = _mm512_roundscale_pd(
                        _mm512_mul_pd(d2dim, ib2),
                        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
                    d2dim = _mm512_fnmadd_pd(t, bx2, d2dim);

                    __m512d d2v = _mm512_mul_pd(d0, d0);
                    d2v = _mm512_fmadd_pd(d1, d1, d2v);
                    d2v = _mm512_fmadd_pd(d2dim, d2dim, d2v);

                    const __m512d r_v =
                        _mm512_mul_pd(_mm512_add_pd(Di_v, Dj), half_v);
                    const __m512d r2_v = _mm512_mul_pd(r_v, r_v);
                    const __mmask8 m = lanemask &
                        _mm512_cmp_pd_mask(d2v, r2_v, _CMP_LT_OQ);
                    tls_candidates += remaining >= 8 ? 8 : remaining;
                    if (swap_edges) {   // SWAP: collect near pairs
                        const __m512d cm_v = _mm512_add_pd(r_v, _mm512_set1_pd(swap_margin));
                        unsigned nm = lanemask & _mm512_cmp_pd_mask(
                            d2v, _mm512_mul_pd(cm_v, cm_v), _CMP_LT_OQ);
                        if (nm) {
                            alignas(64) double d2s[8];
                            _mm512_store_pd(d2s, d2v);
                            while (nm) {
                                const int ln = __builtin_ctz(nm);
                                nm &= nm - 1;
                                swap_edges[swap_n++] = {static_cast<std::uint32_t>(i),
                                                        static_cast<std::uint32_t>(pi[jdx + ln]),
                                                        d2s[ln]};
                            }
                        }
                    }
                    if (!m) continue;

                    alignas(64) double d2_arr[8], r_arr[8];
                    alignas(64) double dx0_arr[8], dx1_arr[8], dx2_arr[8];
                    _mm512_store_pd(d2_arr, d2v);
                    _mm512_store_pd(r_arr, r_v);
                    _mm512_store_pd(dx0_arr, d0);
                    _mm512_store_pd(dx1_arr, d1);
                    _mm512_store_pd(dx2_arr, d2dim);
                    unsigned mm = m;
                    while (mm) {
                        const int lane = __builtin_ctz(mm);
                        mm &= mm - 1;
                        const std::size_t j = pi[jdx + lane];
                        const double r_ij = r_arr[lane];
                        const double d2 = d2_arr[lane];
                        const double dxl[3] = {dx0_arr[lane], dx1_arr[lane],
                                               dx2_arr[lane]};
                        ++tls_overlaps;
                        double dist = std::sqrt(d2);
                        z_t[i]++;
                        z_t[j]++;
                        if (i < tls_lo) tls_lo = i;
                        if (j < tls_lo) tls_lo = j;
                        if (i + 1 > tls_hi) tls_hi = i + 1;
                        if (j + 1 > tls_hi) tls_hi = j + 1;
                        double F_mag = -K * (r_ij / dist - 1.0);
                        double overlap = 1.0 - dist / r_ij;
                        U_red += overlap;
                        ++count_red;
                        Lc_red += r_ij;
                        if (overlap > max_min_red) { max_min_red = overlap; }
                        for (std::size_t d = 0; d < 3; ++d) {
                            double fcomp = F_mag * dxl[d];
                            F_t[i * 3 + d] += fcomp;
                            F_t[j * 3 + d] -= fcomp;
                        }
                        dkappa_red = dkappa_red + K * (rk_data ? (rk_data[i] + rk_data[j]) : r_ij) * (dist - r_ij);
                        Fmean_red = Fmean_red + F_mag;
                        if (swap_g) { swap_g[i] += r_ij - dist; swap_g[j] += r_ij - dist; }   // SWAP proxy
                    }
                }
                // Cycle 21b: masked block covers the tail — nothing left
                // for the scalar loop below.
                jdx = numNbr + 1;
            } else if (use_simd2) {
                // Cycle 21f: 2D variant of the batched reject test — same
                // op-for-op mirroring of the scalar path, two coordinate
                // gathers instead of three.
                const __m512d xi0 = _mm512_set1_pd(xi[0]);
                const __m512d xi1 = _mm512_set1_pd(xi[1]);
                const __m512d Di_v = _mm512_set1_pd(Di);
                const __m512d half_v = _mm512_set1_pd(0.5);
                const __m512d kap_v = _mm512_set1_pd(kappa);
                const double* pk_i =
                    use_packed ? packed_d0 + pair_offsets[i] : nullptr;
                const __m512d ib0 = _mm512_set1_pd(inv_box[0]);
                const __m512d ib1 = _mm512_set1_pd(inv_box[1]);
                const __m512d bx0 = _mm512_set1_pd(box_data[0]);
                const __m512d bx1 = _mm512_set1_pd(box_data[1]);
                const __m256i one2 = _mm256_set1_epi32(2);
                for (std::uint32_t blk = 1; blk <= numNbr; blk += 8) {
                    jdx = blk;
                    const std::uint32_t remaining = numNbr - blk + 1;
                    const __mmask8 lanemask =
                        remaining >= 8
                            ? static_cast<__mmask8>(0xFF)
                            : static_cast<__mmask8>((1u << remaining) - 1u);
                    const __m256i jv = _mm256_maskz_loadu_epi32(
                        lanemask, pi + blk);
                    const __m256i j2 = _mm256_mullo_epi32(jv, one2);
                    const __m512d zero_v = _mm512_setzero_pd();
                    const __m512d xj0 = _mm512_mask_i32gather_pd(
                        zero_v, lanemask, j2, x_data, 8);
                    const __m512d xj1 = _mm512_mask_i32gather_pd(
                        zero_v, lanemask,
                        _mm256_add_epi32(j2, _mm256_set1_epi32(1)), x_data, 8);
                    // Cycle 21k: packed-D0 contiguous load (see 3D variant).
                    const __m512d Dj = use_packed
                        ? _mm512_mul_pd(
                              _mm512_maskz_loadu_pd(lanemask, pk_i + blk),
                              kap_v)
                        : _mm512_mask_i32gather_pd(
                              zero_v, lanemask, jv, D_data, 8);

                    __m512d d0 = _mm512_sub_pd(xj0, xi0);
                    __m512d t = _mm512_roundscale_pd(
                        _mm512_mul_pd(d0, ib0),
                        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
                    d0 = _mm512_fnmadd_pd(t, bx0, d0);

                    __m512d d1 = _mm512_sub_pd(xj1, xi1);
                    t = _mm512_roundscale_pd(
                        _mm512_mul_pd(d1, ib1),
                        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
                    d1 = _mm512_fnmadd_pd(t, bx1, d1);

                    __m512d d2v = _mm512_mul_pd(d0, d0);
                    d2v = _mm512_fmadd_pd(d1, d1, d2v);

                    const __m512d r_v =
                        _mm512_mul_pd(_mm512_add_pd(Di_v, Dj), half_v);
                    const __m512d r2_v = _mm512_mul_pd(r_v, r_v);
                    const __mmask8 m = lanemask &
                        _mm512_cmp_pd_mask(d2v, r2_v, _CMP_LT_OQ);
                    tls_candidates += remaining >= 8 ? 8 : remaining;
                    if (swap_edges) {   // SWAP: collect near pairs
                        const __m512d cm_v = _mm512_add_pd(r_v, _mm512_set1_pd(swap_margin));
                        unsigned nm = lanemask & _mm512_cmp_pd_mask(
                            d2v, _mm512_mul_pd(cm_v, cm_v), _CMP_LT_OQ);
                        if (nm) {
                            alignas(64) double d2s[8];
                            _mm512_store_pd(d2s, d2v);
                            while (nm) {
                                const int ln = __builtin_ctz(nm);
                                nm &= nm - 1;
                                swap_edges[swap_n++] = {static_cast<std::uint32_t>(i),
                                                        static_cast<std::uint32_t>(pi[jdx + ln]),
                                                        d2s[ln]};
                            }
                        }
                    }
                    if (!m) continue;

                    alignas(64) double d2_arr[8], r_arr[8];
                    alignas(64) double dx0_arr[8], dx1_arr[8];
                    _mm512_store_pd(d2_arr, d2v);
                    _mm512_store_pd(r_arr, r_v);
                    _mm512_store_pd(dx0_arr, d0);
                    _mm512_store_pd(dx1_arr, d1);
                    unsigned mm = m;
                    while (mm) {
                        const int lane = __builtin_ctz(mm);
                        mm &= mm - 1;
                        const std::size_t j = pi[jdx + lane];
                        const double r_ij = r_arr[lane];
                        const double d2 = d2_arr[lane];
                        const double dxl[2] = {dx0_arr[lane], dx1_arr[lane]};
                        ++tls_overlaps;
                        double dist = std::sqrt(d2);
                        z_t[i]++;
                        z_t[j]++;
                        if (i < tls_lo) tls_lo = i;
                        if (j < tls_lo) tls_lo = j;
                        if (i + 1 > tls_hi) tls_hi = i + 1;
                        if (j + 1 > tls_hi) tls_hi = j + 1;
                        double F_mag = -K * (r_ij / dist - 1.0);
                        double overlap = 1.0 - dist / r_ij;
                        U_red += overlap;
                        ++count_red;
                        Lc_red += r_ij;
                        if (overlap > max_min_red) { max_min_red = overlap; }
                        for (std::size_t d = 0; d < 2; ++d) {
                            double fcomp = F_mag * dxl[d];
                            F_t[i * 2 + d] += fcomp;
                            F_t[j * 2 + d] -= fcomp;
                        }
                        dkappa_red = dkappa_red + K * (rk_data ? (rk_data[i] + rk_data[j]) : r_ij) * (dist - r_ij);
                        Fmean_red = Fmean_red + F_mag;
                        if (swap_g) { swap_g[i] += r_ij - dist; swap_g[j] += r_ij - dist; }   // SWAP proxy
                    }
                }
                jdx = numNbr + 1;
            }
#endif
            for (; jdx <= numNbr; ++jdx) {
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
                if (swap_edges) {   // SWAP: collect near pairs, branch-free append
                    const double cm = r_ij + swap_margin;
                    swap_edges[swap_n] = {static_cast<std::uint32_t>(i),
                                          static_cast<std::uint32_t>(j), d2};
                    swap_n += (d2 < cm * cm) ? 1u : 0u;
                }

                const double r_ij_sq = r_ij * r_ij;
                if (d2 < r_ij_sq) {
                    ++tls_overlaps;
                    double dist = std::sqrt(d2);
                    z_t[i]++;
                    z_t[j]++;
                    if (i < tls_lo) tls_lo = i;
                    if (j < tls_lo) tls_lo = j;
                    if (i + 1 > tls_hi) tls_hi = i + 1;
                    if (j + 1 > tls_hi) tls_hi = j + 1;
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

                    dkappa_red = dkappa_red + K * (rk_data ? (rk_data[i] + rk_data[j]) : r_ij) * (dist - r_ij);
                    Fmean_red = Fmean_red + F_mag;
                    if (swap_g) { swap_g[i] += r_ij - dist; swap_g[j] += r_ij - dist; }   // SWAP proxy
                }
            }
        }

        if (swap_edges) swap_sink->count[tid] = swap_n;   // SWAP: publish edge count
        if (swap_g) { swap_sink->g_lo[tid] = tls_lo; swap_sink->g_hi[tid] = tls_hi; }  // SWAP: g dirty range
        // Cycle 21e: publish this thread's touched range for the
        // range-limited zero (next call) and reduce (below).
        if (tid < 64) {
            g_fl_lo[tid] = tls_lo;
            g_fl_hi[tid] = tls_hi;
        }
        // Cycle 14 H53: always accumulate candidates/overlaps so we can
        // compute force.loop rejection rate without RCP_CYCLE_PROBE.
        g_floop_candidates.fetch_add(tls_candidates, std::memory_order_relaxed);
        g_floop_overlaps.fetch_add(tls_overlaps, std::memory_order_relaxed);
        if (probe) {
            tls_cycles = rcp_read_cycles() - tsc_start;
            g_floop_cycles.fetch_add(tls_cycles, std::memory_order_relaxed);
            g_floop_test_cycles.fetch_add(tls_cyc_test, std::memory_order_relaxed);
            g_floop_overlap_cycles.fetch_add(tls_cyc_overlap, std::memory_order_relaxed);
        }
    }
    g_t_forces_loop.end();

    // Cycle 21e: record range metadata so the next call's zeroing can be
    // range-limited (done after the loop so all threads have published).
    g_fl_threads = nthreads;
    g_fl_N = N;

    g_t_forces_reduce.begin();
    // Cycle 14: parallel-over-i reduce. Each i sums across all threads. F[i]
    // is written by exactly one thread (the one whose chunk owns i). Order
    // of summation across t = 0..nthreads-1 is fixed for determinism.
    // Cycle 21e: skip threads whose touched range excludes i — their slot
    // holds exactly 0.0 there, and adding 0.0 is the identity (the only
    // representable difference would be a -0.0, which propagates to
    // identical downstream values; byte A/B adjudicates).
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        double Fi_local[8] = {0};
        std::size_t zi = 0;
        for (int tt = 0; tt < nthreads; ++tt) {
            if (i < g_fl_lo[tt] || i >= g_fl_hi[tt]) continue;
            const double*      F_t = F_data + tt * F_stride + i * Ndim;
            const std::size_t* z_t = z_data + tt * z_stride;
            for (std::size_t d = 0; d < Ndim; ++d) Fi_local[d] += F_t[d];
            zi += z_t[i];
        }
        const std::size_t base = i * Ndim;
        for (std::size_t d = 0; d < Ndim; ++d) F[base + d] += Fi_local[d];
        z[i] += zi;
    }
    if (swap_sink && swap_sink->want_gsum) {   // SWAP proxy: reduce per-thread g over
        double* g_out = swap_sink->gsum.data();  // each thread's touched range only
        const std::size_t* glo = swap_sink->g_lo.data();
        const std::size_t* ghi = swap_sink->g_hi.data();
        const double* gbase = swap_sink->gsum_thread.data();
        const std::size_t gstride = swap_sink->gstride;
        #pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < N; ++i) {
            double gi = 0.0;
            for (int tt = 0; tt < nthreads; ++tt) {
                if (i < glo[tt] || i >= ghi[tt]) continue;
                gi += gbase[static_cast<std::size_t>(tt) * gstride + i];
            }
            g_out[i] = gi;
        }
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
                            dkappa = dkappa + 2 * K * (rk_data ? rk_data[i] : r_ij) * (dist - r_ij);
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
                    dkappa = dkappa + 2 * K * (rk_data ? rk_data[i] : r_ij) * (d_ij - r_ij);
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

    dkappa = dkappa + mu * (sumD_active >= 0.0 ? sumD_active : sumD);  // cascade: active-only pressure
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

// ---- ellipsoids: glue between the public structs and rcpgenerator/ellipsoid.hpp (additive) -------------------
namespace {

// The rule at the core input. Returns the shape and fills the flat rotation array (N * Ndim * Ndim).
// Is this the non-spherical path: aspect ratios (ellipsoids), rod inputs, or orientations?
bool ell_input_present(const PackingInput& input, const PackingConfig& config) {
    return !config.aspect_ratios.empty() || !input.aspect_ratios.empty() || !input.orientations.empty() ||
           !input.rod_length.empty() || !input.rod_radius.empty();
}

// Force pass of the non-spherical path: ellipsoids or rods, by the shape's kind.
void shape_forces(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x, std::size_t N,
                  const EllShape& sh, const std::vector<double>& Db, const std::vector<double>& box,
                  const std::vector<std::int8_t>& walls, const double* Rot, double* F, double* Tq, double& U,
                  double& max_min_dist, double& Lc, double& Fmean, double mu, double& dkappa,
                  std::vector<std::size_t>& z, double& energy, std::vector<double>& Amat_work,
                  EllPairCache* cache = nullptr, double* gsum = nullptr, RodPairCache* rod_cache = nullptr) {
    if (sh.kind != 0) {
        rod_forces(pairs_data, pair_offsets, x, N, sh, Db, box, walls, Rot, F, Tq, U, max_min_dist, Lc, Fmean, mu, dkappa,
                   z, energy, Amat_work, gsum, rod_cache);
        return;
    }
    ell_forces(pairs_data, pair_offsets, x, N, sh, Db, box, walls, Rot, F, Tq, U, max_min_dist, Lc, Fmean, mu, dkappa, z,
               energy, Amat_work, cache, gsum);
}

void ell_read_orientations(const PackingInput& input, std::size_t N, std::size_t Ndim, std::vector<double>& rot,
                           const char* who) {
    if (input.orientations.size() != N) {
        throw std::runtime_error(std::string(who) + ": orientations needs one rotation per particle (" +
                                 std::to_string(N) + "), got " + std::to_string(input.orientations.size()));
    }
    rot.assign(N * Ndim * Ndim, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        if (input.orientations[i].size() != Ndim * Ndim) {
            throw std::runtime_error(std::string(who) + ": each orientation is a row-major Ndim x Ndim rotation matrix (" +
                                     std::to_string(Ndim * Ndim) + " numbers); particle " + std::to_string(i) +
                                     " has " + std::to_string(input.orientations[i].size()));
        }
        for (std::size_t k = 0; k < Ndim * Ndim; ++k) rot[i * Ndim * Ndim + k] = input.orientations[i][k];
        ell_check_orientation(rot.data() + i * Ndim * Ndim, Ndim, i);
    }
}

EllShape ell_resolve_input(const PackingInput& input, const PackingConfig& config, std::size_t N,
                           std::size_t Ndim, std::vector<double>& rot, std::vector<double>& qtab,
                           std::vector<double>& deq, std::vector<double>& bfac,
                           std::vector<double>* rod_user = nullptr) {
    bfac.assign(N, 1.0);   // bounding diameter / user's D: 1 for ellipsoids (D is the longest axis)
    const bool rods = !input.rod_length.empty() || !input.rod_radius.empty();
    if (rods) {
        if (!config.aspect_ratios.empty() || !input.aspect_ratios.empty()) {
            throw std::runtime_error("rods: rod_length / rod_radius were supplied together with aspect_ratios; a packing "
                                     "is of rods or of ellipsoids, not both");
        }
        if (input.orientations.empty()) {
            throw std::runtime_error("rods: rod_length and rod_radius were supplied without orientations; supply both "
                                     "(rods) or neither (spheres). At the Packing level the initializer generates "
                                     "orientations for you");
        }
        EllShape sh;
        sh.Ndim = Ndim;
        sh.kind = config.rod_capped != 0 ? SHAPE_ROD_CAPPED : SHAPE_ROD_UNCAPPED;
        sh.n_iter = std::max(config.ell_iters, 30);   // bound of the 3D uncapped iteration (most pairs need none)
        // sh.track_tol stays 0 here (single evaluations iterate to convergence); a run sets it from ell_track_tol
        sh.all_one = false;
        sh.rod_round = config.rod_round > 0.0 ? config.rod_round : 0.0;   // single evaluations; a quench run sets its own
        rod_build_shape_table(input.rod_length, input.rod_radius, config.rod_capped != 0, N, Ndim, qtab, deq, bfac);
        if (rod_user) {   // the user's own numbers per particle (length, radius), handed back unchanged at the end
            rod_user->assign(2 * N, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                (*rod_user)[2 * i] = input.rod_length.size() == 1 ? input.rod_length[0] : input.rod_length[i];
                (*rod_user)[2 * i + 1] = input.rod_radius.size() == 1 ? input.rod_radius[0] : input.rod_radius[i];
            }
        }
        sh.qtab = qtab.data();
        sh.deqtab = deq.data();
        ell_read_orientations(input, N, Ndim, rot, "rods");
        return sh;
    }
    if (config.aspect_ratios.empty() && input.aspect_ratios.empty()) {
        throw std::runtime_error("ellipsoids: orientations were supplied without aspect_ratios (or rod_length and "
                                 "rod_radius); supply both (ellipsoids or rods) or neither (spheres)");
    }
    if (input.orientations.empty()) {
        throw std::runtime_error("ellipsoids: aspect_ratios were supplied without orientations; supply both "
                                 "(ellipsoids) or neither (spheres). At the Packing level the initializer "
                                 "generates orientations for you");
    }
    // One list for every particle or a table of N rows; either way the kernel works from the per-particle table
    // (so a table of equal rows and the single list are the same computation).
    EllShape sh;
    sh.Ndim = Ndim;
    sh.n_iter = config.ell_iters < 1 ? 1 : config.ell_iters;
    ell_build_shape_table(config.aspect_ratios, input.aspect_ratios, N, Ndim, qtab, deq, sh.all_one);
    sh.qtab = qtab.data();
    sh.deqtab = deq.data();
    if (input.orientations.size() != N) {
        throw std::runtime_error("ellipsoids: orientations needs one rotation per particle (" +
                                 std::to_string(N) + "), got " + std::to_string(input.orientations.size()));
    }
    rot.assign(N * Ndim * Ndim, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
        if (input.orientations[i].size() != Ndim * Ndim) {
            throw std::runtime_error("ellipsoids: each orientation is a row-major Ndim x Ndim rotation matrix (" +
                                     std::to_string(Ndim * Ndim) + " numbers); particle " + std::to_string(i) +
                                     " has " + std::to_string(input.orientations[i].size()));
        }
        for (std::size_t k = 0; k < Ndim * Ndim; ++k) rot[i * Ndim * Ndim + k] = input.orientations[i][k];
        ell_check_orientation(rot.data() + i * Ndim * Ndim, Ndim, i);
    }
    return sh;
}

std::vector<std::vector<double>> ell_nested_orientations(const std::vector<double>& rot, std::size_t N,
                                                         std::size_t Ndim) {
    std::vector<std::vector<double>> out(N, std::vector<double>(Ndim * Ndim, 0.0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t k = 0; k < Ndim * Ndim; ++k) out[i][k] = rot[i * Ndim * Ndim + k];
    return out;
}

// Mean |torque| over particles carrying one, relative to (mean fractional contact force x mean contact length):
// the rotational counterpart of the reported force magnitude. 0 for spheres.
double ell_torque_magnitude(const std::vector<double>& tq, std::size_t N, std::size_t Ndim, double Fmean,
                            double Lc) {
    const std::size_t nrot = ell_nrot(Ndim);
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < N; ++i) {
        double t2 = 0.0;
        for (std::size_t c = 0; c < nrot; ++c) t2 += tq[i * nrot + c] * tq[i * nrot + c];
        if (t2 > 1.0e-32) { sum += std::sqrt(t2); ++count; }
    }
    if (count == 0 || !(Fmean > 0.0) || !(Lc > 0.0)) return 0.0;
    return (sum / double(count)) / (Fmean * Lc * Lc);
}

// All-pairs half list in the layout the force kernels read: row i = [count, j > i ...].
void all_pairs_layout(std::size_t N, std::vector<std::uint64_t>& offsets, std::vector<std::uint32_t>& data) {
    offsets.assign(N + 1, 0);
    for (std::size_t i = 0; i < N; ++i) offsets[i + 1] = offsets[i] + (N - i);
    data.assign(offsets[N] + 1, 0);
    for (std::size_t i = 0; i < N; ++i) {
        std::uint32_t* row = data.data() + offsets[i];
        row[0] = static_cast<std::uint32_t>(N - 1 - i);
        for (std::size_t j = i + 1; j < N; ++j) row[j - i] = static_cast<std::uint32_t>(j);
    }
}

struct FlatInput {
    std::size_t N = 0, Ndim = 0;
    std::vector<double> x, D, box;
    std::vector<std::int8_t> walls;
};

FlatInput flatten_input(const PackingInput& input, const PackingConfig& config) {
    FlatInput f;
    f.N = input.positions.size();
    f.Ndim = f.N ? input.positions[0].size() : 0;
    if (f.N == 0 || f.Ndim < 2) throw std::runtime_error("need at least one particle in Ndim >= 2");
    if (input.diameters.size() != f.N) throw std::runtime_error("diameters length does not match positions");
    f.x.assign(f.N * f.Ndim, 0.0);
    for (std::size_t i = 0; i < f.N; ++i)
        for (std::size_t d = 0; d < f.Ndim; ++d) f.x[i * f.Ndim + d] = input.positions[i][d];
    f.D = input.diameters;
    f.box = config.box.empty() ? std::vector<double>(f.Ndim, 1.0) : config.box;
    f.walls = config.walls.empty() ? std::vector<std::int8_t>(f.Ndim, 0) : config.walls;
    if (f.box.size() != f.Ndim || f.walls.size() != f.Ndim) throw std::runtime_error("box/walls length mismatch");
    return f;
}

}  // namespace

SingleEvaluation evaluate_once(const PackingInput& input, const PackingConfig& config, double mu) {
    FlatInput f = flatten_input(input, config);
    const std::size_t N = f.N, Ndim = f.Ndim;
    const bool ell = ell_input_present(input, config);
    std::vector<std::uint64_t> offsets;
    std::vector<std::uint32_t> pairs;
    all_pairs_layout(N, offsets, pairs);
    std::vector<double> F(N * Ndim, 0.0);
    SingleEvaluation out;
    std::size_t nrot = ell_nrot(Ndim);
    out.torques.assign(N, std::vector<double>(nrot, 0.0));
    if (ell) {
        std::vector<double> rot, amat, tq(N * nrot, 0.0), qtab, deq, bfac;
        const EllShape sh = ell_resolve_input(input, config, N, Ndim, rot, qtab, deq, bfac);
        for (std::size_t i = 0; i < N; ++i) f.D[i] *= bfac[i];
        shape_forces(pairs.data(), offsets.data(), f.x.data(), N, sh, f.D, f.box, f.walls, rot.data(), F.data(),
                   tq.data(), out.U, out.max_overlap, out.Lc, out.Fmean, mu, out.dkappa, out.z, out.energy, amat);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t c = 0; c < nrot; ++c) out.torques[i][c] = tq[i * nrot + c];
    } else {
        std::vector<std::vector<double>> min_dist;
        get_forces_nd_3(pairs.data(), offsets.data(), f.x.data(), N, Ndim, f.D, f.box, f.walls, F.data(), out.U,
                        min_dist, out.max_overlap, out.Lc, out.Fmean, mu, out.dkappa, out.z);
        // The sphere kernel never evaluates an energy; the one its forces derive from is written out here.
        double e = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                double d2 = 0.0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double delta = f.x[j * Ndim + d] - f.x[i * Ndim + d];
                    if (f.walls[d] == 0) delta -= std::nearbyint(delta / f.box[d]) * f.box[d];
                    d2 += delta * delta;
                }
                const double r_ij = 0.5 * (f.D[i] + f.D[j]);
                if (d2 < r_ij * r_ij) { const double ov = r_ij - std::sqrt(d2); e += 0.5 * ov * ov; }
            }
            for (std::size_t d = 0; d < Ndim; ++d) {
                if (f.walls[d] != 1) continue;
                for (int wall = 0; wall < 2; ++wall) {
                    const double dist = std::abs(f.box[d] * wall - f.x[i * Ndim + d]);
                    if (dist < 0.5 * f.D[i]) e += (0.5 * f.D[i] - dist) * (0.5 * f.D[i] - dist);
                }
            }
            if (f.walls[0] < 0) {   // curved container: force 2 (r - gap) towards the axis, energy (r - gap)^2
                const std::size_t tc = curved_axes(f.walls[0], Ndim);
                const double Rc = 0.5 * f.box[0];
                double dc2 = 0.0;
                for (std::size_t m = 0; m < tc; ++m) { const double u = f.x[i * Ndim + m] - Rc; dc2 += u * u; }
                const double gap = Rc - std::sqrt(dc2);
                if (gap < 0.5 * f.D[i]) e += (0.5 * f.D[i] - gap) * (0.5 * f.D[i] - gap);
            }
        }
        out.energy = e;
    }
    out.forces.assign(N, std::vector<double>(Ndim, 0.0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t d = 0; d < Ndim; ++d) out.forces[i][d] = F[i * Ndim + d];
    return out;
}

ContactScale contact_scale(const PackingInput& input, const PackingConfig& config) {
    FlatInput f = flatten_input(input, config);
    const std::size_t N = f.N, Ndim = f.Ndim;
    const bool ell = ell_input_present(input, config);
    EllShape sh;
    std::vector<double> rot, qtab, deq, bfac;
    if (ell) {
        sh = ell_resolve_input(input, config, N, Ndim, rot, qtab, deq, bfac);
        for (std::size_t i = 0; i < N; ++i) f.D[i] *= bfac[i];
        sh.vol_factor = sh.kind != 0 ? rod_volume_factor(N, Ndim, f.D.data(), qtab.data(), sh.kind == SHAPE_ROD_CAPPED)
                                     : ell_volume_factor(N, Ndim, f.D.data(), qtab.data());
    } else {
        sh = make_ell_shape(std::vector<double>(Ndim - 1, 1.0), Ndim, config.ell_iters);
        rot.assign(N * Ndim * Ndim, 0.0);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t d = 0; d < Ndim; ++d) rot[i * Ndim * Ndim + d * Ndim + d] = 1.0;
    }
    ContactScale out;
    out.scale = sh.kind != 0
        ? rod_min_contact_scale(N, sh, f.x.data(), f.D.data(), rot.data(), f.box, f.walls, out.pair_i, out.pair_j, out.n_overlapping)
        : ell_min_contact_scale(N, sh, f.x.data(), f.D.data(), rot.data(), f.box, f.walls, out.pair_i,
                                      out.pair_j, out.n_overlapping);
    double vol = 0.0;
    for (double dv : f.D) vol += sphere_volume(0.5 * dv, Ndim);
    out.particle_volume = vol * sh.vol_factor;
    return out;
}

namespace {
PermuteOnce ell_permute_outside(const PackingInput& input, const PackingConfig& config, bool forced,
                                std::size_t fi, std::size_t fj) {
    FlatInput f = flatten_input(input, config);
    const std::size_t N = f.N, Ndim = f.Ndim, nn = Ndim * Ndim;
    std::vector<double> rot, qtab, deq, bfac, rodu;
    EllShape sh = ell_resolve_input(input, config, N, Ndim, rot, qtab, deq, bfac, &rodu);
    for (std::size_t i = 0; i < N; ++i) f.D[i] *= bfac[i];
    auto swap_user = [&](std::size_t a, std::size_t b) { if (!rodu.empty()) { std::swap(rodu[2 * a], rodu[2 * b]); std::swap(rodu[2 * a + 1], rodu[2 * b + 1]); } };
    PermuteOnce out;
    std::vector<EllExchangeRecord> done;
    const bool move = config.ell_permute_orient != 0;
    if (forced) {
        if (fi >= N || fj >= N || fi == fj) throw std::runtime_error("exchange_once: i and j must be two different particles");
        ell_exchange(Ndim, fi, fj, f.D.data(), nullptr, qtab.data(), deq.data(), rot.data(), move);
        std::swap(bfac[fi], bfac[fj]); swap_user(fi, fj);
    } else {
        std::vector<std::uint64_t> offsets;
        std::vector<std::uint32_t> pairs;
        all_pairs_layout(N, offsets, pairs);
        std::vector<double> amat(N * nn);
        if (sh.kind == 0) ell_shape_matrices(N, sh, f.D.data(), rot.data(), amat.data());
        EllPermuteParams ep;
        ep.window = config.swap_geo_w;
        ep.budget = config.swap_geo_B > 0 ? config.swap_geo_B : static_cast<long>(N);
        ep.max_partners = config.ell_permute_partners > 0 ? config.ell_permute_partners : static_cast<int>(N);
        ep.move_orientations = move;
        ep.decide_in_place = config.ell_permute_orient != 1;
        ep.overlap_guard = config.ell_permute_guard != 0;
        ep.n_iter = sh.n_iter;
        ep.kind = sh.kind;
        out.candidates = ell_permute_pass(pairs.data(), offsets.data(), f.x.data(), N, Ndim, f.box, f.walls,
                                          f.D.data(), nullptr, qtab.data(), deq.data(), rot.data(), amat.data(), ep, done);
        for (const EllExchangeRecord& e : done) { std::swap(bfac[e.i], bfac[e.j]); swap_user(e.i, e.j); }
    }
    for (std::size_t i = 0; i < N; ++i) f.D[i] /= bfac[i];
    if (sh.kind != 0) {   // rods: the shape out as length and radius per particle
        out.rod_length.assign(N, 0.0); out.rod_radius.assign(N, 0.0);
        for (std::size_t i = 0; i < N; ++i) { out.rod_length[i] = rodu[2 * i]; out.rod_radius[i] = rodu[2 * i + 1]; }
    }
    out.diameters = f.D;
    out.orientations = ell_nested_orientations(rot, N, Ndim);
    out.aspect_ratios.assign(N, std::vector<double>(Ndim - 1, 1.0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t k = 1; k < Ndim; ++k) out.aspect_ratios[i][k - 1] = qtab[i * Ndim + k];
    for (const EllExchangeRecord& e : done) {
        out.exchanges.push_back(static_cast<double>(e.i)); out.exchanges.push_back(static_cast<double>(e.j)); out.exchanges.push_back(e.gain);
    }
    return out;
}
}  // namespace

PermuteOnce permute_once(const PackingInput& input, const PackingConfig& config) {
    return ell_permute_outside(input, config, false, 0, 0);
}
PermuteOnce exchange_once(const PackingInput& input, const PackingConfig& config, std::size_t i, std::size_t j) {
    return ell_permute_outside(input, config, true, i, j);
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
    // ---- ellipsoids (additive; rcpgenerator/ellipsoid.hpp) --------------------------------------------------
    // The rule at the core input: no aspect ratios and no orientations -> ell_on == false and nothing below
    // that is guarded by ell_on runs (the sphere path, bit for bit); both -> the ellipsoid path; one without
    // the other -> error. In the ellipsoid path D holds the BOUNDING diameter (D_user x rmax) from here to
    // result assembly, so the neighbour machinery works unchanged; see the header for the conventions.
    const bool ell_on = ell_input_present(input, config);
    EllShape ell_shape;
    std::vector<double> ell_rot, ell_torque, ell_m, ell_v, ell_amat;
    std::vector<double> ell_rodu;         // rods: the user's (length, radius) per particle, to hand back bit for bit
    std::vector<double> ell_bfac;         // bounding diameter / user's D per particle (1 for ellipsoids; rods: rod.hpp)
    std::vector<double> ell_q, ell_deq;   // per-particle shape (N x Ndim) and volume-equivalent factor; follow the particle
    auto ell_sync = [&]() { ell_shape.qtab = ell_q.data(); ell_shape.deqtab = ell_deq.data(); };
    std::vector<EllExchangeRecord> ell_exch;   // permutation: the exchanges of the current pass
    std::size_t ell_perm_candidates = 0;
    long ell_perm_productive = 0;              // passes that made at least one exchange
    bool ell_capture_armed = false;
    long ell_cap_step = -1; double ell_cap_energy = 0.0;
    std::vector<double> ell_cap_x, ell_cap_D, ell_cap_rot, ell_cap_q, ell_cap_F, ell_cap_T, ell_cap_exch, ell_cap_z;
    double ell_q_stall_u = -1.0;  // quench: energy at the stall rule's last reference point
    std::vector<double> ell_qw;   // quench: rotational FIRE velocities, in length units (angular velocity x bounding radius)
    EllPairCache ell_cache;   // each pair's lam from the previous step
    RodPairCache rod_cache;   // 3D uncapped rods: each pair's contact normal from the previous step
    double ell_energy = 0.0;
    if (ell_on) {
        ell_shape = ell_resolve_input(input, config, N, Ndim, ell_rot, ell_q, ell_deq, ell_bfac, &ell_rodu);
        ell_shape.track_tol = config.ell_track_tol > 0.0 ? config.ell_track_tol : 0.0;
        if (config.quench_mode != 0) {
            // A quench balances forces to ~1e-10 of the contact load, so the contact function is converged
            // (lam iterated until its step is below 1e-15, started from the pair's previous value) instead of
            // tracked to ell_track_tol: the forces are then the exact gradient of the energy, to rounding.
            ell_shape.track_tol = 1e-300;
            ell_shape.conv_tol = 1e-15;
            ell_shape.n_iter = std::max(ell_shape.n_iter, 30);
        }
        ell_cache.reset(N);
        if (ell_shape.kind != 0) ell_shape.rod_round = config.quench_mode != 0 && config.quench_rod_round > 0.0 ? config.quench_rod_round : 0.0;
        if (ell_shape.kind == SHAPE_ROD_UNCAPPED && Ndim == 3 && !(ell_shape.rod_round > 0.0)) rod_cache.reset(N);
        for (std::size_t i = 0; i < N; ++i) D[i] *= ell_bfac[i];
        const std::size_t nrot = ell_nrot(Ndim);
        ell_torque.assign(N * nrot, 0.0);
        ell_m.assign(N * nrot, 0.0);
        ell_v.assign(N * nrot, 0.0);
        if (config.quench_mode != 0) ell_qw.assign(N * nrot, 0.0);
    }
    // Cycle 14: x and x_last are flat [N*Ndim]. Convert from the caller's
    // nested input once here; convert back to nested at result time.
    std::vector<double> x(N * Ndim, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t d = 0; d < Ndim; ++d)
            x[i * Ndim + d] = input.positions[i][d];
    std::vector<double> x_last = x;
    // Cascade (additive): per-particle frozen-diameter mask. Empty/all-zero ->
    // has_fixed=false and every cascade branch below is skipped, so the one-shot
    // path is byte-identical. Permuted alongside D0 on reorder (see below).
    std::vector<std::int8_t> fixed_flags = input.fixed;
    if (!fixed_flags.empty() && fixed_flags.size() != N) fixed_flags.assign(N, 0);
    const bool has_fixed = std::any_of(fixed_flags.begin(), fixed_flags.end(),
                                       [](std::int8_t f) { return f != 0; });
    // Cascade (additive): the frozen-diameter mask is a STARTUP scaffold only.
    // At the rung-1->0 trigger it is released (see below) so every class is free
    // to densify / shrink-to-contact together (McGeary final settle). Guarded by
    // has_fixed throughout -> the one-shot path (has_fixed=false) never sees it.
    bool cascade_released = false;
    // Cycle 20: per-particle diameter anchor for the growth-aware refresh
    // trigger. Diameter growth between refreshes consumes pair-shell budget
    // exactly like displacement does ((Di+Dj)/2 grows while the candidate
    // list stays fixed); the sentinel caught real missed pairs from this in
    // the 2D early phase where kappa grows fastest. Every site that resets
    // x_last must reset D_last alongside it.
    std::vector<double> D_last = D;
    // Cycle 20/21h: refresh-trigger scale, resolved once per run.
    //   theta = scale * Dmin/4.
    // Safety: a pair can close by up to 3*theta between refreshes (theta
    // from the freshly-rebuilt side, 2*theta from the side whose anchor
    // predates the rebuild), so 3*theta <= t_shell is required, where the
    // pair-finders use t_shell = min(max(median*1.02, Dmin*2.25), Dmax)/2.
    // Cycle 21h: the scale is now DERIVED from that bound for any
    // distribution (no special cases): with rho = t_shell/Dmin
    // (scale-invariant under the uniform kappa rescale, so computed once),
    //   theta = 0.9 * rho*Dmin/3  ->  scale = 1.2*rho   (10% margin).
    // For the broad distributions used in production rho = 1.125 -> scale
    // 1.35 (unchanged behavior); for monodisperse rho = 0.5 -> scale 0.6.
    // The generality suite's missed-pair sentinel caught the monodisperse
    // bound violation in the ORIGINAL fixed Dmin/4 trigger (t_shell
    // collapses to Dmax/2 < 0.75*Dmin there) — this formula closes it.
    // Env override RCP_REFRESH_TRIGGER_SCALE (capped 2.0) takes precedence;
    // values above the bound are only safe with a matching
    // RCP_T_SHELL_SCALE raise.
    const double trigger_scale_run = [&]() {
        const char* s = std::getenv("RCP_REFRESH_TRIGGER_SCALE");
        if (s && s[0]) {
            const double v = std::atof(s);
            if (v > 0.0 && v <= 2.0) return v;
        }
        double dmx = 0.0;
        double dmn = std::numeric_limits<double>::max();
        for (double dv : D) {
            if (dv > dmx) dmx = dv;
            if (dv < dmn) dmn = dv;
        }
        std::vector<double> Ds = D;
        std::nth_element(Ds.begin(), Ds.begin() + Ds.size() / 2, Ds.end());
        const double med = Ds[Ds.size() / 2];
        const double t_shell =
            std::min(std::max(med * 1.02, dmn * 2.25), dmx) * 0.5;
        const double rho = t_shell / dmn;
        // Cap at 1.9: keeps theta = scale*Dmin/4 below the kdtree
        // patch-path bbox inflation floor (Dmin/2) with margin, for
        // mid-heavy distributions where 1.2*rho would exceed it. The cap
        // is strictly more conservative than the 3*theta bound.
        // Cycle 21k: 1.2 -> 1.28 (margin 10% -> 4%). The zero-margin bound
        // limit is (4/3)*rho; 1.29*rho was sentinel-clean in the Cycle 21
        // sweep (so was the zero-margin 1.33*rho) and cut refreshed
        // particles ~10% / pairs time ~8%. The running-min growth
        // accounting (21j) makes consumed budget exact, which is what
        // justifies thinning the margin.
        return std::min(1.28 * rho, 1.9);
    }();

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
    if (ell_on) {
        // Features not yet extended to ellipsoids refuse aspect ratios (they must never run on the sphere
        // assumption).
        const bool rods_on = ell_shape.kind != 0;
        auto refuse = [rods_on](const char* what) {
            if (rods_on) throw std::runtime_error(std::string("rods (rod_length / rod_radius) are not supported with ") + what +
                                                  " yet; remove the rod inputs or do not use that feature");
            throw std::runtime_error(std::string("ellipsoids (aspect_ratios) are not supported with ") + what +
                                     " yet; remove aspect_ratios or do not use that feature");
        };
        // Diameter permutation is extended to ellipsoids in the form the public controls use (capacity-matched
        // exchanges on a fixed cadence: swap_geo == 1 with the proxy sweeps off). The engine's research variants
        // work from sphere quantities and are not.
        if (config.swap > 0.0 &&
            (config.swap_geo != 1 || !config.swap_sweeps_off || config.swap_exact != 0 || config.swap_rescue > 0.0 ||
             config.swap_cycle3 > 0.0 || config.swap_geo_T > 0.0 || config.swap_geo_cycle4 != 0 ||
             config.swap_geo_capk != 0 || config.swap_geo_cage != 0 || config.swap_geo_bottleneck != 0 ||
             config.swap_geo_nonsep != 0 || config.swap_geo_reopen != 0 || config.swap_event_log_stride > 0)) {
            refuse("the permutation engine's research variants (only the diameter_permute recipe is extended)");
        }
        // The quench is extended to ellipsoids (positions and orientations); its research-only variants are not.
        if (config.quench_mode != 0 && (config.quench_two_level != 0 || config.quench_freeze_damp > 0.0 ||
                                        config.quench_warm_steps > 0 || config.quench_mix_mode > 1)) {
            refuse("the quench's two-level, activity-damping, warm-up or per-particle-reset variants");
        }
        if (!config.morph_target.empty()) refuse("the distribution morph");
        if (has_fixed) refuse("fixed (frozen-diameter) particles");
        if (!config.swap_inject.empty()) refuse("swap injection");
    }
    // Ellipsoids: every volume in the loop is the TRUE particle volume = bounding-sphere volume x ell_vf, so
    // phi, phi_history and every threshold on them mean the same thing as for spheres. ell_vf is used only
    // inside ell_on branches; the sphere path never multiplies by it.
    // With a shape per particle ell_vf is the ratio of the two totals (true volume / bounding-sphere volume):
    // one number for the run, because rescaling, reordering and exchanging whole particles leave it unchanged.
    const double ell_vf = !ell_on ? 1.0
        : ell_shape.kind != 0 ? rod_volume_factor(N, Ndim, D.data(), ell_q.data(), ell_shape.kind == SHAPE_ROD_CAPPED)
        : ell_volume_factor(N, Ndim, D.data(), ell_q.data());

    // [226] fix_height: the box arrives REALIZED (the initializer already set box[Ndim-1] = multiple x the largest
    // diameter, and both Packing and run_packing pass that box on). Multiplying by D[0] again here shrank the height
    // below one particle diameter (and D[0] is the largest particle only before the first spatial reorder).

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
        auto scaled_pair = scale_diameters_nd(D, ell_on ? phi * phi_modifier / ell_vf : phi * phi_modifier, box, Ndim, fix_height);
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
    if (ell_on) current_volume *= ell_vf;

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
        // [227] Ndim >= 4 with the AUTOMATIC base (neighbor_max = 0): the table base of 5500 reserves
        // ~7,400 (4D) to ~14,800 (8D) slots per particle against a measured need of 39 / 79 / 150
        // (max, monodisperse jammed, 4D / 5D / 6D) — at N = 2M in 5D that is ~74 GB of address space,
        // more than the allocator grants on typical machines. Instead the cap is 4 x the expected
        // number of neighbours inside the search shell of a jammed packing, plus 100: centres inside
        // the reach ball R_i = (D_i + 1.01 Dmax)/2 + t, i.e. phi_J(d) R_i^d / <(D/2)^d>, with
        // phi_J(d) = 0.448 x 0.66^(d-4) (fit to the measured 0.448 / 0.300 / 0.193). Scale-free and
        // per particle (large particles get more). A user-supplied neighbor_max keeps the rule above;
        // 2D and 3D are unchanged. Capacity never changes results, only whether a run fits.
        double est_coef = 0.0, est_dmax = 0.0, est_t = 0.0;
        if (Ndim >= 4 && N > 0 && config.neighbor_max == 0) {
            est_dmax = *std::max_element(D.begin(), D.end());
            std::vector<double> Ds(D.begin(), D.end());
            std::nth_element(Ds.begin(), Ds.begin() + N / 2, Ds.end());
            est_t = std::min(std::max(Ds[N / 2] * 1.02, dmin * 2.25), est_dmax) * 0.5;
            double mean_vol = 0.0;
            for (std::size_t i = 0; i < N; ++i) mean_vol += std::pow(D[i] * 0.5, static_cast<double>(Ndim));
            mean_vol /= static_cast<double>(N);
            const double phi_j = 0.448 * std::pow(0.66, static_cast<double>(Ndim) - 4.0);
            if (mean_vol > 0.0) est_coef = 4.0 * phi_j / mean_vol;
        }
        max_nb_out.assign(N, 0);
        offsets_out.assign(N + 1, 0);
        std::uint64_t off = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const double ratio = D[i] / dmin;
            double cap_d = scaling *
                (50.0 * std::pow(ratio, omega) + static_cast<double>(base));
            if (est_coef > 0.0) {
                const double reach = (D[i] + 1.01 * est_dmax) * 0.5 + est_t;
                cap_d = est_coef * std::pow(reach, static_cast<double>(Ndim)) + 100.0;
                cap_d = std::min(cap_d, static_cast<double>(N));   // before the integer cast
            }
            std::uint32_t cap = static_cast<std::uint32_t>(std::ceil(cap_d));
            if (static_cast<std::uint64_t>(cap) > N) {
                cap = static_cast<std::uint32_t>(N);
            }
            max_nb_out[i] = cap;
            offsets_out[i] = off;
            off += static_cast<std::uint64_t>(cap) + 1;
        }
        offsets_out[N] = off;
        if (std::getenv("RCP_DEBUG_LAYOUT")) {
            std::uint32_t cmin = max_nb_out[0], cmax = max_nb_out[0];
            for (std::size_t i = 1; i < N; ++i) {
                if (max_nb_out[i] < cmin) cmin = max_nb_out[i];
                if (max_nb_out[i] > cmax) cmax = max_nb_out[i];
            }
            std::cerr << "[layout] N=" << N << " base=" << base
                      << " dmin=" << dmin << " caps=[" << cmin << ","
                      << cmax << "] total_slots=" << off << "\n";
        }
        return off;
    };

    std::vector<std::uint32_t> max_neighbors_per_particle;
    std::vector<std::uint64_t> pair_offsets;
    rcp_rss_probe("pre-pair-layout");
    std::uint64_t pairs_total = compute_pair_layout(
        max_neighbors_per_particle, pair_offsets);
    // Cycle 21L lazy pairs buffer. The N x M layout is CAPACITY — for
    // mid-heavy size distributions occupancy is a few percent (e.g.
    // gaussian bulk ~20x above a tail dmin: ~5000 slots/particle
    // allocated, ~60 used). A zero-filled vector touches every page and
    // makes the full 2+ GB capacity RESIDENT for the whole run; with
    // many search workers running fat candidates simultaneously that
    // was a ~100 GB host-memory wave. Untouched anonymous pages are
    // zero-backed by the kernel, and row content beyond slot[0]=count
    // is never read — so allocate WITHOUT value-init and reset only the
    // per-row count slots. RSS then tracks actual usage, not capacity.
    struct LazyPairs {
        std::unique_ptr<std::uint32_t[]> buf;
        std::uint64_t n = 0;
        std::uint32_t* data() { return buf.get(); }
        std::uint64_t size() const { return n; }
        std::uint32_t* begin() { return buf.get(); }
        std::uint32_t* end() { return buf.get() + n; }
        std::uint32_t& operator[](std::uint64_t i) { return buf[i]; }
        const std::uint32_t& operator[](std::uint64_t i) const { return buf[i]; }
        void reset(std::uint64_t total, const std::uint64_t* offsets,
                   std::size_t N) {
            if (total != n) {
                buf.reset(new std::uint32_t[total]);   // no value-init
                n = total;
            }
            for (std::size_t i = 0; i < N; ++i)
                buf[offsets[i]] = 0;                   // counts only
        }
        // Cycle 23: return faulted pages to the OS. SAFE ONLY when the whole
        // buffer is about to be fully rewritten (the reorder-triggered full
        // rebuild): old contents are discarded anyway, and any page touched
        // again zero-fills then gets overwritten. Without this, each reorder
        // recomputes pair_offsets (per-particle caps shift even though the
        // total is invariant), so the rebuild writes a SHIFTED layout into the
        // same buffer and faults fresh pages while the prior layout's pages
        // stay resident — RSS ratchets to the full-matrix footprint over the
        // run. Only acts on the large mmap-backed (page-aligned) allocation.
        void release_resident() {
#if defined(__linux__)
            const std::size_t bytes = n * sizeof(std::uint32_t);
            if (!buf || bytes < (std::size_t{64} << 20)) return;  // large mmap-backed only
            constexpr std::uintptr_t PAGE = 4096;
            // malloc returns the user ptr offset into a page-aligned mmap region
            // (chunk header), so it is 16-aligned, NOT page-aligned. Round the
            // start UP to the next page and madvise whole pages of the interior.
            const std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(buf.get());
            const std::uintptr_t aligned = (addr + PAGE - 1) & ~(PAGE - 1);
            const std::size_t off = static_cast<std::size_t>(aligned - addr);
            if (bytes <= off) return;
            const std::size_t len = ((bytes - off) / PAGE) * PAGE;
            if (!len) return;
            const int rc = ::madvise(reinterpret_cast<void*>(aligned), len,
                                     MADV_DONTNEED);
            static const bool dbg = []() {
                const char* s = std::getenv("RCP_DEBUG_LAYOUT");
                return s && s[0] && s[0] != '0';
            }();
            if (dbg)
                std::fprintf(stderr, "[release_resident] dropped %.2f GB rc=%d\n",
                             len / 1e9, rc);
#endif
        }
    };
    LazyPairs pairs_data;
    pairs_data.reset(pairs_total, pair_offsets.data(), N);
    rcp_rss_probe("post-pairs-alloc");
    // Cycle 21k packed-D0 mirror: per-slot D0[j] copies, maintained by the
    // pair backends' canonicalize passes, consumed by the force loop's
    // SIMD reject paths (contiguous load + kappa multiply instead of a
    // D[j] gather; byte-identical by the single-rounding argument). The
    // helper self-heals after any pairs_data resize — every resize site
    // is followed by a full-refresh dispatch, which rebuilds all rows.
    std::vector<double> packed_d0_mirror;
    // Cascade: the packed-D0 SIMD mirror multiplies a per-slot D0 by the global
    // kappa, which would wrongly scale frozen particles. Disable it when has_fixed
    // so the force loop falls back to the per-particle D[j] gather (correct sizes).
    const bool packed_d0_on = rcp_packed_d0() && (Ndim == 2 || Ndim == 3) && !has_fixed;
    auto packed_d0_ptr = [&]() -> double* {
        if (!packed_d0_on) return nullptr;
        if (packed_d0_mirror.size() != pairs_data.size())
            packed_d0_mirror.assign(pairs_data.size(), 0.0);
        return packed_d0_mirror.data();
    };
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
    //
    // Cycle 22+: termination is governed by convergence, not by a fixed step
    // ceiling — annealing mu -> 0 drives the maximum overlap below the criterion,
    // which is guaranteed regardless of step count (a low-mu state simply cannot
    // retain overlaps). The old flat 60000 default cap was cutting off the
    // floor-mu relaxation for the largest-N cases, leaving them under-converged.
    // We replace it with a schedule-derived ceiling that gives provable headroom
    // over a full run: (rung-1 deadline) x (oscillation cycles + 3), floored at
    // N_STEPS for backward compatibility. It stays finite so the step-indexed
    // history buffers below remain bounded (an unbounded cap would re-trigger the
    // OOB corruption noted above), and options.max_steps > 0 is still honored.
    //
    // NOTE: the preset table below MIRRORS the in-loop schedule (RCP_SPEED ->
    // (Y, X, cycles) and the rung-1 deadline, ~lines 5180+ and ~5755). Keep the
    // two in sync if the schedule constants change.
    const std::size_t step_cap = [&]() -> std::size_t {
        if (config.max_steps > 0)                                   // fail-soft cap [log 043]
            return static_cast<std::size_t>(config.max_steps);
        if (options.max_steps > 0) return options.max_steps;        // passable override
        int Y = 2, X = 9; std::size_t cycles = 2;                   // RCP_SPEED=quick default
        if (const char* sp = std::getenv("RCP_SPEED")) {
            std::string v(sp);
            if      (v == "immediate") { Y = 1; X = 5;  cycles = 1; }
            else if (v == "patient")   { Y = 3; X = 24; cycles = 3; }
            else if (v == "forever")   { Y = 4; X = 49; cycles = 4; }
        }
        if (const char* oc = std::getenv("RCP_OSC_CYCLES")) {
            long v = std::atol(oc); if (v >= 1 && v < 100) cycles = std::size_t(v);
        }
        double rung1 = static_cast<double>(Y) * 4000.0
            + static_cast<double>(X) * std::sqrt(static_cast<double>(N)
                                                 * static_cast<double>(Ndim))
            - 2000.0;                                               // Cycle 21M -2000 cut
        if (const char* rd = std::getenv("RCP_RUNG1_DEADLINE")) {
            long v = std::atol(rd);
            if (v > 100 && v < 10000000) rung1 = static_cast<double>(v);
        }
        if (rung1 < 1500.0) rung1 = 1500.0;
        const std::size_t formula =
            static_cast<std::size_t>(rung1 * static_cast<double>(cycles + 3));
        return std::max<std::size_t>(static_cast<std::size_t>(N_steps), formula);
    }();
    const std::size_t hist_size = step_cap + 1;
    std::vector<double> U_history(hist_size, 0.0);

    std::vector<double> phi_history(hist_size, 0.0);
    std::vector<double> F_history(hist_size, 0.0);
    // Cycle 16: ADAM convergence study — capture max overlap, mu, alpha,
    // mu_flag per step for offline diagnostic. Same size as phi/F histories.
    std::vector<double> max_overlap_history(hist_size, 0.0);
    std::vector<double> mu_hist(hist_size, 0.0);
    std::vector<double> alpha_hist(hist_size, 0.0);
    std::vector<std::int8_t> mu_flag_hist(hist_size, 0);
    std::vector<float> zpos_hist(config.record_visibility ? hist_size : 0, 0.f);
    std::vector<float> mueff_hist(config.record_visibility ? hist_size : 0, 0.f);

    double phi_min = 0.8;
    if (Ndim == 2) { phi_min = 0.76; }
    if (Ndim == 3) { phi_min = 0.45; }
    if (Ndim == 4) { phi_min = 0.15; }
    if (Ndim == 5) { phi_min = 0.1; }
    if (Ndim > 5) { phi_min = 0.1 / std::pow(2, (Ndim - 4)); }

    double dkappa = 1.0;
    double kappa = 1.0;
    // Cycle 21k adaptive shell: per-particle trigger budget, anchored at
    // each particle's refresh with the shell scale its list was actually
    // built under (theta_k = Dmin/4 * trigger_scale * s_at_build). The
    // bound 3*theta <= t_shell holds per-list because the backends widen
    // t_shell by the same s_at_build. With RCP_SHELL_SMAX=1 every anchor
    // stores exactly the old dmin_quarter and behavior is byte-identical.
    std::vector<double> theta_anchor(N, Dmin * 0.25 * trigger_scale_run);
    double shell_gr_ema = 0.0;
    double shell_kappa_prev = 1.0;
    // Bulk re-anchor: every full-refresh site (reorder, state load,
    // target-phi lock, Dmin-jump rebuild, rung -1 force-refresh) resets
    // x_last/D_last for ALL particles; budgets re-anchor with them.
    auto shell_refill_anchor = [&]() {
        theta_anchor.assign(N,
            Dmin * 0.25 * trigger_scale_run * g_rcp_shell_scale_dyn);
    };
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
    const double mu1_ref = mu;                 // mu-hold reference (rung-1 mu) [log 071]
    long mu_hold_applied = 0;
    double mu_ctrl_mult = config.mu_hold_mult;  // visibility controller state [log 073]
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
        if (ell_on) s *= ell_vf;   // ellipsoids: S_D0 carries the true volume
        return std::make_pair(dmin, s);
    };
    auto d0_cache = recompute_d0_cache();
    double D0_min_cached = d0_cache.first;
    double S_D0_cached   = d0_cache.second;
    // Cycle 21k: D0 max for the adaptive-shell controller (worst-case
    // per-step budget burn is on the largest particle).
    auto recompute_d0_max = [&]() {
        double m = 0.0;
        for (double v : D0) if (v > m) m = v;
        return m;
    };
    double D0_max_cached = recompute_d0_max();

    // Cycle 10 identity tracking: particle_origin[k] tells us which original
    // index the particle currently at slot k started at. Updated alongside
    // every reorder. At the end of the run, used to compare against an
    // un-reordered baseline particle-by-particle.
    std::vector<std::size_t> particle_origin(N);
    for (std::size_t k = 0; k < N; ++k) particle_origin[k] = k;

    // Diameter-SWAP (additive). Everything guarded by swap_on; when config.swap == 0
    // no swap code touches any array the packing loop reads (constraint C4).
    const bool swap_on = config.swap > 0.0;
    SwapParams swap_params;
    swap_params.swap = config.swap;
    swap_params.interval = std::max(1, config.swap_interval);
    swap_params.mode = config.swap_mode;
    swap_params.rank_window = config.swap_rank_window;
    swap_params.size_tol = config.swap_size_tol;
    swap_params.eps = config.swap_eps;
    swap_params.start = config.swap_start;
    swap_params.stop = config.swap_stop;
    swap_params.phase = config.swap_phase;
    swap_params.adam = config.swap_adam;
    swap_params.max_per_particle = config.swap_max_per_particle;
    swap_params.exact = config.swap_exact;
    swap_params.window_rel = config.swap_window_rel;
    swap_params.rank_space = config.swap_rank_space;
    swap_params.ctrl = config.swap_ctrl;
    swap_params.target_a = config.swap_a;
    swap_params.target_R = config.swap_R;
    swap_params.ctrl_eta = config.swap_eta;
    swap_params.w_min = config.swap_wmin;
    swap_params.w_max = config.swap_wmax;
    swap_params.s_max = config.swap_smax;
    swap_params.floor_budget = config.swap_floor_budget;
    swap_params.rescue = config.swap_rescue;
    swap_params.rescue_w = config.swap_rescue_w;
    swap_params.rescue_interval = config.swap_rescue_interval;
    swap_params.cycle3 = config.swap_cycle3;
    swap_params.event_log_stride = config.swap_event_log_stride;
    swap_params.event_log_cap = config.swap_event_log_cap;
    swap_params.event_log_geom = config.swap_event_log_geom;
    swap_params.geo = config.swap_geo;
    swap_params.geo_T = config.swap_geo_T;
    swap_params.geo_T_end = config.swap_geo_T_end;       // [135] H9 annealing
    swap_params.geo_T_decay = config.swap_geo_T_decay;
    swap_params.geo_w = config.swap_geo_w;
    swap_params.geo_B = config.swap_geo_B;
    swap_params.geo_k = std::max<long>(1, config.swap_geo_k);
    swap_params.geo_start = config.swap_geo_start;
    swap_params.geo_stop = config.swap_geo_stop;
    swap_params.geo_cycle = config.swap_geo_cycle;             // H2 [117]
    swap_params.geo_cycle_cap = config.swap_geo_cycle_cap;
    swap_params.geo_cycle4 = config.swap_geo_cycle4;           // H34 [194]
    swap_params.geo_cycle4_cap = config.swap_geo_cycle4_cap;
    swap_params.geo_capk = config.swap_geo_capk;               // H7 [121/122]
    swap_params.geo_cage = config.swap_geo_cage;               // H8 [129/130/131]
    swap_params.geo_bottleneck = config.swap_geo_bottleneck;   // H12 [146/148]
    swap_params.geo_cage_move = config.swap_geo_cage_move;     // [133] diagnostic position move
    swap_params.geo_nonsep = config.swap_geo_nonsep;          // [172] H19 non-separable comparator
    swap_params.geo_ctarget = config.swap_geo_ctarget;
    swap_params.geo_nsk = config.swap_geo_nsk;
    swap_params.geo_nonsep_start = config.swap_geo_nonsep_start;   // [190] H26 window
    swap_params.geo_nonsep_stop = config.swap_geo_nonsep_stop;
    swap_params.geo_capmode = config.swap_geo_capmode;
    swap_params.idle_frac = config.swap_idle_frac;
    SwapEngine swap_engine(swap_params, static_cast<std::uint64_t>(config.seed));
    std::vector<std::uint32_t> swap_touched;
    std::vector<double> swap_log_buf;      // per sweep: [step, mu_flag, cand, acc, dE, phi, contacts, e_scale]
    std::uint64_t swap_total_accepted = 0;
    SwapEdgeSink swap_sink;                // near-pair edges collected by get_forces_nd_3 on sweep steps
    std::vector<double> swap_D0_layout;    // D0 at the last pair-layout computation (capacity guard)
    bool swap_need_relayout = false;
    std::size_t swap_relayouts = 0;
    std::size_t swap_next_step = static_cast<std::size_t>(std::max<long>(0, swap_params.start));
    std::size_t swap_cur_interval = static_cast<std::size_t>(swap_params.interval);
    // [219] DISTRIBUTION MORPH state (inert unless morph_target is set): rank-keyed quantile interpolation
    // of D0 from its state at morph_start toward the (sorted, volume-matched) target; feedback on max overlap.
    const bool morph_on = !config.morph_target.empty();
    const bool quench_on = config.quench_mode != 0;                      // [221]
    const bool mu_pinned = config.mu_pin > 0.0 || quench_on;
    std::vector<double> q_v(quench_on ? N * Ndim : 0, 0.0);              // FIRE velocities
    std::vector<char> q_held(N, 0);
    std::vector<char> q_light(N, 0);    // [223] two-level: 1 = relaxed by the local overdamped step, excluded from FIRE's sums
    double q_dt = config.quench_dt0, q_alpha = config.quench_alpha0; long q_npos = 0, q_iter = 0;
    std::vector<double> q_act; double q_act_med = 0.0;   // [224] activity-gated dissipation state
    double q_res = 1.0, q_gres = 1.0; bool q_converged = false;
    // [225] give-up triggers and wall bookkeeping. Wall contacts enter F in the force kernel but not z or the
    // per-particle contact-force sums, so the quench adds them here; with no walls none of this runs.
    int q_status = 0; long q_unjam_streak = 0, q_stall_ref_iter = 0; double q_stall_ref = -1.0; bool q_any_held = false;
    std::vector<double> q_wall_s; std::vector<int> q_wall_n;
    std::vector<double> q_log; std::vector<double> q_held_force(Ndim, 0.0); std::vector<double> q_minv;
    bool morph_armed = false;               // reference captured at the first update step
    double morph_lambda = 0.0, morph_rate = config.morph_rate0, morph_ref_maxov = 0.0, morph_vol0 = 0.0;
    long morph_ok_streak = 0, morph_updates = 0, morph_holds = 0;
    double morph_rate_ceiling = config.morph_rate_max;   // lowered to the last hold-free sustained rate on every trip
    double morph_kappa0 = 0.0;                            // tracking baseline for the deflation guard (EMA over morph_kappa_tau updates)
    double morph_last_safe = 0.0;                         // largest increment sustained for morph_grow_after in-band updates
    long morph_hold_streak = 0;                           // consecutive held updates (kappa-lag stall detector)
    bool morph_stalled = false;
    std::vector<double> morph_q_start, morph_q_target;   // ascending quantile tables (size N)
    std::vector<std::uint32_t> morph_order;              // argsort of D0 (ascending), rebuilt when the permutation epoch changes
    std::uint64_t morph_order_epoch = ~0ULL;              // epoch at which morph_order was built (reorders/trades bump it)
    std::vector<double> morph_log_buf;                   // 5 per update: step, lambda, increment, maxov, held

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
        // Cycle 21d: deterministic D-stat cache invalidation — the kdtree's
        // cached min/max/median indices are meaningless after a permutation.
        ++g_rcp_perm_epoch;
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
        // Cascade: the frozen-diameter mask must travel with the spatial reorder.
        if (!fixed_flags.empty()) {
            std::vector<std::int8_t> tmp(N);
            for (std::size_t k = 0; k < N; ++k) tmp[k] = fixed_flags[perm[k]];
            fixed_flags = std::move(tmp);
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
        if (ell_on) {   // ellipsoids: orientation and its Adam moments follow their particle
            auto apply_perm_stride = [&](std::vector<double>& arr, std::size_t stride) {
                std::vector<double> tmp(N * stride);
                for (std::size_t k = 0; k < N; ++k)
                    for (std::size_t d = 0; d < stride; ++d) tmp[k * stride + d] = arr[perm[k] * stride + d];
                arr = std::move(tmp);
            };
            apply_perm_stride(ell_rot, Ndim * Ndim);
            apply_perm_stride(ell_q, Ndim);
            apply_perm_stride(ell_deq, 1);
            apply_perm_stride(ell_bfac, 1);
            if (!ell_rodu.empty()) apply_perm_stride(ell_rodu, 2);
            ell_sync();
            apply_perm_stride(ell_m, ell_nrot(Ndim));
            apply_perm_stride(ell_v, ell_nrot(Ndim));
            if (!ell_qw.empty()) apply_perm_stride(ell_qw, ell_nrot(Ndim));
            ell_cache.clear();   // indices changed: every pair starts cold on the next force pass
            if (!rod_cache.count.empty()) rod_cache.clear();
        }
        if (!q_v.empty()) apply_perm_flat(q_v);          // [221] FIRE velocities follow their particles through the reorder
        if (quench_on) {
            std::vector<char> th(N); for (std::size_t k = 0; k < N; ++k) th[k] = q_held[perm[k]]; q_held.swap(th);
            { std::vector<char> tl(N); for (std::size_t k = 0; k < N; ++k) tl[k] = q_light[perm[k]]; q_light.swap(tl); }
            if (q_minv.size() == N) { std::vector<double> tm(N); for (std::size_t k = 0; k < N; ++k) tm[k] = q_minv[perm[k]]; q_minv.swap(tm); }
        }
        // Cycle 17 H6 fix: reorder triggers full pair-list rebuild via
        // `refresh = [1,1,...]` at line below, so K-batched accumulated
        // refresh flags are no longer relevant — zero them. Original code
        // left them pointing at old slot indices, causing pair-list churn
        // and F-spike on the FIRST post-reorder get_pairs call.
        std::fill(accumulated_refresh.begin(), accumulated_refresh.end(), 0u);
        accumulated_has_refresh = false;

        // Cycle 13: D has been permuted, so per-particle MaxNb caps and
        // offsets shift too. Recompute the layout and resize pairs_data.
        // Cycle 23: drop the prior layout's resident pages first — the full
        // rebuild below (refresh=[1,...]) overwrites the whole buffer, so this
        // is result-neutral and stops RSS ratcheting across shifting layouts.
        pairs_data.release_resident();
        std::uint64_t new_total = compute_pair_layout(
            max_neighbors_per_particle, pair_offsets);
        pairs_data.reset(new_total, pair_offsets.data(), N);

        // Rebuild pair list fresh under the new numbering.
        std::fill(refresh.begin(), refresh.end(), 1u);
        g_t_pairs.begin();
        {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min,
                               packed_d0_ptr(), D0.data());
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
        g_t_pairs.end();
        x_last = x;
        D_last = D;
        shell_refill_anchor();   // Cycle 21k
        // SWAP: index<->particle mapping changed; rank tables are rebuilt lazily and
        // the per-particle swap counters travel with the particles. The layout was
        // just recomputed from the current D -> snapshot D0 for the capacity guard.
        if (swap_on) {
            swap_D0_layout = D0;
            swap_engine.invalidate_ranks();
            auto& sc = swap_engine.swap_counts();
            if (sc.size() == N) {
                std::vector<std::uint64_t> tmp(N);
                for (std::size_t k = 0; k < N; ++k) tmp[k] = sc[perm[k]];
                sc = std::move(tmp);
            }
        }
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
    g_t_swap.reset();
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
    // Cycle 21d: new run = new particle indexing — invalidate the kdtree
    // D-stat cache (it is static and would otherwise survive across runs
    // in the same process).
    ++g_rcp_perm_epoch;
    // Cycle 21k: the adaptive shell scale is static — reset per run.
    g_rcp_shell_scale_dyn = 1.0;
    // checkpoint-deterministic mode for this run [log 044/B4]
    g_rcp_ckpt_deterministic = !config.checkpoint_path.empty();
    auto profile_wall_start = std::chrono::steady_clock::now();

    rcp_rss_probe("pre-init-reorder");
    // Cycle 10: initial spatial reorder so the very first force loop sees
    // a cache-friendly layout. Skipped if REORDER_INTERVAL=0 (baseline).
    if (options.reorder_interval > 0) {
        reorder_particles_now();
        rcp_rss_probe("post-init-reorder");
    } else {
        g_t_pairs.begin();
        {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min,
                               packed_d0_ptr(), D0.data());
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
        g_t_pairs.end();
    }

    g_t_forces.begin();
    if (ell_on) {
        shape_forces(pairs_data.data(), pair_offsets.data(), x.data(), N, ell_shape, D, box, walls, ell_rot.data(),
                   F, ell_torque.data(), U, max_min_dist, Lc, Fmean, mu, dkappa, z, ell_energy, ell_amat,
                   ell_shape.track_tol > 0.0 ? &ell_cache : nullptr, nullptr, rod_cache.count.empty() ? nullptr : &rod_cache);
    } else
    get_forces_nd_3(pairs_data.data(), pair_offsets.data(), x.data(), N, Ndim, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);
    g_t_forces.end();

    PackingTrace trace;
    std::size_t traj_log_next = 1, traj_last_captured = 0;   // log-spaced samples: next step to take; last step sampled
    auto record_observation = [&](std::size_t step_value, bool capture_sample_positions) {
        PackingObservation observation;
        observation.step = step_value;
        observation.phi = phi / phi_modifier;
        observation.force_magnitude = F_magnitude;
        observation.energy = U;
        observation.max_min_dist = max_min_dist;

        if (capture_sample_positions) {
            trace.steps.push_back(step_value);
            if (ell_on) {   // ellipsoids: report the user's D (first axis) and the orientations
                std::vector<double> Du(D);
                for (std::size_t i = 0; i < N; ++i) Du[i] /= ell_bfac[i];
                trace.diameters.push_back(std::move(Du));
                trace.orientations.push_back(ell_nested_orientations(ell_rot, N, Ndim));
            } else
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
                trace.origins.push_back(particle_origin);
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
    // Number of steps whose per-step history entries are fully written (all
    // *_history series filled through [0, hist_filled)). The history buffers are
    // over-allocated to hist_size = step_cap + 1 for indexing headroom; this
    // tracks the actually-filled prefix so the returned series carry no trailing
    // zero (see the resize before the result is assembled).
    std::size_t hist_filled = 0;
    const std::size_t max_steps = step_cap;   // resolved above: N-scaled default (cap no longer binds real runs), options.max_steps honored if passed
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
    // The force threshold of the stop test. Flat-ended rods have their own (Ken, 2026-11-13; FINDINGS C39): their
    // normalised force does not fall below about 6e-3 (2D) and 1e-2 to 5e-2 (3D) in this scheme, whatever the
    // tightness of the contact solve (flat contacts; LOG [038]), so 5e-3 could never be met. 2e-2 and 1e-1 lie
    // above the largest force seen in the 750 steps before the stop in 18 + 18 survey cases: the force test
    // rarely decides when these shapes stop. An explicit RCP_OSC_F_THRESH applies to every shape.
    const double stop_f_thresh =
        (ell_on && ell_shape.kind == SHAPE_ROD_UNCAPPED && !(std::getenv("RCP_OSC_F_THRESH") && std::getenv("RCP_OSC_F_THRESH")[0]))
            ? (Ndim == 2 ? 2e-2 : 1e-1) : osc_f_thresh;
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
            D_last = D;
            shell_refill_anchor();   // Cycle 21k
            D0 = ls.D0;
            m_flat = ls.m_flat;
            v_flat = ls.v_flat;
            // d0_cache, D0_min_cached, S_D0_cached depend on D0; recompute.
            d0_cache = recompute_d0_cache();
            D0_min_cached = d0_cache.first;
            S_D0_cached = d0_cache.second;
            D0_max_cached = recompute_d0_max();   // Cycle 21k
            // Pair-layout caps (MaxNb_i and pair_offsets) are computed from
            // D's actual range; the initial call at line ~3637 used the
            // pre-load D which has a different scale. Recompute now so the
            // first get_pairs call after load has correct slot allocation.
            {
                std::uint64_t new_total = compute_pair_layout(
                    max_neighbors_per_particle, pair_offsets);
                pairs_data.reset(new_total, pair_offsets.data(), N);
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
    if (swap_on && swap_D0_layout.size() != N) swap_D0_layout = D0;   // capacity-guard baseline

    // ===================== checkpoint v3 [log 043] ==================================
    // Periodic bitwise state checkpoints + resume. Opt-in via config.checkpoint_path.
    // Layout: header (magic RCPS, ver 3, N, Ndim, step, seed, step_cap), a fixed
    // 16-double inspector block (analysis/ckpt_inspect.py mirrors it), then the full
    // opaque state in the exact order below (save/load are symmetric; a mismatch is
    // a hard error, never a silent skip). Pair rows are dumped SPARSELY (count +
    // occupied slots) so the lazy capacity buffer's untouched pages stay untouched.
    const bool ckpt_on = !config.checkpoint_path.empty();
    if (ckpt_on && (capture_positions || trajectory_interval > 0))
        throw std::runtime_error("[ckpt] checkpointing does not support trajectory capture");
    std::size_t ckpt_resumed_from = 0;
    std::size_t ckpt_last_step_written = 0;
    auto ckpt_wall_last = std::chrono::steady_clock::now();
    auto ckpt_write = [&](std::size_t at_step) -> bool {
        if (!ckpt_on) return true;
        const std::string tmpp = config.checkpoint_path + ".tmp";
        std::FILE* fp = std::fopen(tmpp.c_str(), "wb");
        if (!fp) { std::cerr << "[ckpt] fopen failed: " << tmpp << "\n"; return false; }
        bool ok = true;
        auto W = [&](const void* q, std::size_t n) { ok = ok && std::fwrite(q, 1, n, fp) == n; };
        auto WV = [&](const auto& vec) {
            const std::uint64_t n = vec.size(); W(&n, 8);
            if (n) W(vec.data(), n * sizeof(vec[0]));
        };
        auto WPFX = [&](const auto& vec, std::size_t used) {
            const std::uint64_t n = std::min<std::uint64_t>(used, vec.size()); W(&n, 8);
            if (n) W(vec.data(), static_cast<std::size_t>(n) * sizeof(vec[0]));
        };
#define RCP_CK(v) W(&(v), sizeof(v))
        if (x.size() != N * Ndim || D0.size() != N || D.size() != N ||
            box.size() != Ndim || particle_origin.size() != N) {
            std::fclose(fp);
            std::cerr << "[ckpt] REFUSING to write: state arrays moved-out/mis-sized"
                      << " (x=" << x.size() << " D=" << D.size() << " D0=" << D0.size()
                      << " box=" << box.size() << " origin=" << particle_origin.size()
                      << ")\n";
            return false;
        }
        const char magic[5] = "RCPS"; const std::uint32_t ver = 3;
        const std::uint64_t N64 = N, Nd64 = Ndim, st64 = at_step,
            seed64 = static_cast<std::uint64_t>(config.seed), cap64 = step_cap;
        W(magic, 4); RCP_CK(ver); RCP_CK(N64); RCP_CK(Nd64); RCP_CK(st64);
        RCP_CK(seed64); RCP_CK(cap64);
        const double insp[16] = { double(at_step), double(mu_flag), mu, alpha, kappa,
            phi / phi_modifier, max_min_dist, F_magnitude, double(int(osc_phase)),
            double(hist_filled), double(swap_total_accepted), double(swap_relayouts),
            shell_gr_ema, F_ema, double(steps), 0.0 };
        W(insp, sizeof(insp));
        RCP_CK(t); RCP_CK(alpha); RCP_CK(mu); RCP_CK(kappa); RCP_CK(m_kappa);
        RCP_CK(v_kappa); RCP_CK(v_update_kappa); RCP_CK(m_hat_kappa); RCP_CK(v_hat_kappa);
        RCP_CK(phi); RCP_CK(coeff); RCP_CK(current_volume); RCP_CK(box_volume);
        RCP_CK(Dmin); RCP_CK(F_magnitude); RCP_CK(max_min_dist); RCP_CK(F_ema);
        RCP_CK(shell_gr_ema); RCP_CK(shell_kappa_prev); RCP_CK(g_rcp_shell_scale_dyn);
        RCP_CK(mu_change); RCP_CK(LastAlphaUpdate); RCP_CK(mu_flag);
        RCP_CK(mu_decay_active); RCP_CK(mu_decay_target); RCP_CK(mu_decay_next_tick_step);
        RCP_CK(alpha_rung1_frozen); RCP_CK(phi_max.first); RCP_CK(phi_max.second);
        RCP_CK(osc_phase); RCP_CK(osc_mu_0); RCP_CK(osc_mu_minus1); RCP_CK(osc_start_step);
        RCP_CK(osc_period); RCP_CK(osc_t_relax); RCP_CK(osc_final_start_step);
        RCP_CK(osc_consec_conv);
        RCP_CK(steps); RCP_CK(hist_filled);
        RCP_CK(swap_next_step); RCP_CK(swap_cur_interval); RCP_CK(swap_need_relayout);
        RCP_CK(swap_relayouts); RCP_CK(swap_total_accepted);
        RCP_CK(cascade_released); RCP_CK(runtime_fix_diameter); RCP_CK(target_phi_locked);
        RCP_CK(accumulated_has_refresh); RCP_CK(update_flag);
        RCP_CK(d0_cache.first); RCP_CK(d0_cache.second);
        RCP_CK(D0_min_cached); RCP_CK(S_D0_cached); RCP_CK(D0_max_cached);
        RCP_CK(agg_max_neighbors); RCP_CK(agg_min_neighbors);
        RCP_CK(g_rcp_perm_epoch); RCP_CK(pp_window_steps);
        RCP_CK(g_dstat_max_idx); RCP_CK(g_dstat_min_idx); RCP_CK(g_dstat_med_idx);
        RCP_CK(g_dstat_N); RCP_CK(g_dstat_D0); RCP_CK(g_dstat_Dmid);
        RCP_CK(g_dstat_valid); RCP_CK(g_dstat_epoch);
        WV(box); WV(x); WV(x_last); WV(D); WV(D0); WV(D_last);
        WV(m_flat); WV(v_flat);
        const std::uint8_t has_aux = rcp_fused_adam() ? 0 : 1; RCP_CK(has_aux);
        if (has_aux) { WV(v_update_flat); WV(m_hat_flat); WV(v_hat_flat); }
        const std::uint8_t has_verlet = (method != "ADAM") ? 1 : 0; RCP_CK(has_verlet);
        if (has_verlet)
            for (auto* nv : {&v_max, &a, &v_verlet, &a_old})
                for (const auto& row : *nv) WV(row);
        WV(fixed_flags); WV(theta_anchor); WV(refresh); WV(accumulated_refresh);
        WV(particle_origin); WV(swap_D0_layout); WV(swap_log_buf); WV(spike_log_buf);
        const std::uint8_t has_pp = pp_enabled ? 1 : 0; RCP_CK(has_pp);
        if (has_pp) { WV(pp_F_sq_sum); WV(pp_F_sign_flips); WV(pp_pos_min_run);
                      WV(pp_pos_max_run); WV(pp_contact_sum); WV(pp_F_sign_prev);
                      WV(pp_signflip_resets); }
        const std::uint8_t has_b6a = b6a_probe_on ? 1 : 0; RCP_CK(has_b6a);
        if (has_b6a) WV(b6a_x_prev_step);
        WPFX(phi_history, hist_filled); WPFX(F_history, hist_filled);
        WPFX(U_history, hist_filled); WPFX(max_overlap_history, hist_filled);
        WPFX(mu_hist, hist_filled); WPFX(alpha_hist, hist_filled);
        WPFX(mu_flag_hist, hist_filled);
        WV(pair_offsets); WV(max_neighbors_per_particle);
        const std::uint64_t pcap = pairs_data.size(); RCP_CK(pcap);
        for (std::size_t i = 0; i < N; ++i) {
            const std::uint32_t cnt = pairs_data[pair_offsets[i]];
            W(&cnt, 4);
            if (cnt) W(&pairs_data[pair_offsets[i] + 1],
                       static_cast<std::size_t>(cnt) * 4);
        }
        const std::uint8_t has_pd0 =
            (packed_d0_on && !packed_d0_mirror.empty()) ? 1 : 0;
        RCP_CK(has_pd0);
        if (has_pd0) WV(packed_d0_mirror);
        if (ell_on) { WV(ell_rot); WV(ell_m); WV(ell_v); WV(ell_q); WV(ell_deq);   // ellipsoids only: sphere blobs are unchanged
                      WV(ell_cache.ids); WV(ell_cache.lam); WV(ell_cache.count); }
        if (ell_on && ell_shape.kind != 0) { WV(ell_bfac); WV(ell_rodu); WV(rod_cache.ids); WV(rod_cache.count); WV(rod_cache.tag); WV(rod_cache.nvec); }   // rods only: ellipsoid blobs are unchanged
        ok = ok && swap_engine.save_state(fp);
        const char tail[5] = "RCPE"; W(tail, 4);
#undef RCP_CK
        std::fclose(fp);
        // debug/archival mode: negative interval keeps every write as a
        // distinct file path.STEP (used by the resume gates to snapshot the
        // unbroken run mid-flight)
        const std::string outp = (config.checkpoint_interval_s < 0 || config.checkpoint_every_steps > 0)
            ? config.checkpoint_path + "." + std::to_string(at_step)
            : config.checkpoint_path;
        if (ok) ok = (std::rename(tmpp.c_str(), outp.c_str()) == 0);
        if (ok) {
            ckpt_last_step_written = at_step;
            if (std::FILE* hb = std::fopen((config.checkpoint_path + ".hb").c_str(), "a")) {
                std::fprintf(hb, "%llu %llu %ld %d %.9g %.9g %.9g %.6g\n",
                             (unsigned long long)at_step,
                             (unsigned long long)step_cap,
                             (long)std::time(nullptr), mu_flag,
                             phi / phi_modifier, F_magnitude, max_min_dist, mu);
                std::fclose(hb);
            }
        } else {
            std::cerr << "[ckpt] WRITE FAILED at step " << at_step << "\n";
        }
        return ok;
    };
    auto ckpt_load = [&]() {
        std::FILE* fp = std::fopen(config.checkpoint_path.c_str(), "rb");
        if (!fp) return;   // no blob yet: fresh run
        bool ok = true;
        auto R = [&](void* q, std::size_t n) { ok = ok && std::fread(q, 1, n, fp) == n; };
        auto RV = [&](auto& vec) {
            std::uint64_t n = 0; R(&n, 8);
            vec.resize(static_cast<std::size_t>(n));
            if (n) R(vec.data(), static_cast<std::size_t>(n) * sizeof(vec[0]));
        };
        auto RPFX = [&](auto& vec) {
            std::uint64_t n = 0; R(&n, 8);
            if (n > vec.size()) { ok = false; return; }
            if (n) R(vec.data(), static_cast<std::size_t>(n) * sizeof(vec[0]));
        };
#define RCP_CK(v) R(&(v), sizeof(v))
        char magic[4]; std::uint32_t ver = 0;
        std::uint64_t N64 = 0, Nd64 = 0, st64 = 0, seed64 = 0, cap64 = 0;
        R(magic, 4); RCP_CK(ver);
        if (!ok || std::memcmp(magic, "RCPS", 4) != 0 || ver != 3)
            throw std::runtime_error("[ckpt] bad magic/version in " + config.checkpoint_path);
        RCP_CK(N64); RCP_CK(Nd64); RCP_CK(st64); RCP_CK(seed64); RCP_CK(cap64);
        if (N64 != N || Nd64 != Ndim ||
            seed64 != static_cast<std::uint64_t>(config.seed))
            throw std::runtime_error("[ckpt] fingerprint mismatch (N/Ndim/seed): "
                                     + config.checkpoint_path);
        // step_cap is a SESSION parameter, not physics: a max_steps-capped run and
        // its uncapped continuation legitimately differ here (fail-soft + resume).
        if (cap64 != step_cap)
            std::cerr << "[ckpt] note: step_cap differs (saved " << cap64
                      << ", now " << step_cap << ") - session parameter, continuing\n";
        double insp[16]; R(insp, sizeof(insp)); (void)insp;
        RCP_CK(t); RCP_CK(alpha); RCP_CK(mu); RCP_CK(kappa); RCP_CK(m_kappa);
        RCP_CK(v_kappa); RCP_CK(v_update_kappa); RCP_CK(m_hat_kappa); RCP_CK(v_hat_kappa);
        RCP_CK(phi); RCP_CK(coeff); RCP_CK(current_volume); RCP_CK(box_volume);
        RCP_CK(Dmin); RCP_CK(F_magnitude); RCP_CK(max_min_dist); RCP_CK(F_ema);
        RCP_CK(shell_gr_ema); RCP_CK(shell_kappa_prev); RCP_CK(g_rcp_shell_scale_dyn);
        RCP_CK(mu_change); RCP_CK(LastAlphaUpdate); RCP_CK(mu_flag);
        RCP_CK(mu_decay_active); RCP_CK(mu_decay_target); RCP_CK(mu_decay_next_tick_step);
        RCP_CK(alpha_rung1_frozen); RCP_CK(phi_max.first); RCP_CK(phi_max.second);
        RCP_CK(osc_phase); RCP_CK(osc_mu_0); RCP_CK(osc_mu_minus1); RCP_CK(osc_start_step);
        RCP_CK(osc_period); RCP_CK(osc_t_relax); RCP_CK(osc_final_start_step);
        RCP_CK(osc_consec_conv);
        RCP_CK(steps); RCP_CK(hist_filled);
        RCP_CK(swap_next_step); RCP_CK(swap_cur_interval); RCP_CK(swap_need_relayout);
        RCP_CK(swap_relayouts); RCP_CK(swap_total_accepted);
        RCP_CK(cascade_released); RCP_CK(runtime_fix_diameter); RCP_CK(target_phi_locked);
        RCP_CK(accumulated_has_refresh); RCP_CK(update_flag);
        RCP_CK(d0_cache.first); RCP_CK(d0_cache.second);
        RCP_CK(D0_min_cached); RCP_CK(S_D0_cached); RCP_CK(D0_max_cached);
        RCP_CK(agg_max_neighbors); RCP_CK(agg_min_neighbors);
        std::uint64_t saved_epoch = 0;
        static_assert(sizeof(g_rcp_perm_epoch) == 8, "epoch width");
        RCP_CK(saved_epoch); RCP_CK(pp_window_steps);
        RCP_CK(g_dstat_max_idx); RCP_CK(g_dstat_min_idx); RCP_CK(g_dstat_med_idx);
        RCP_CK(g_dstat_N); RCP_CK(g_dstat_D0); RCP_CK(g_dstat_Dmid);
        RCP_CK(g_dstat_valid); RCP_CK(g_dstat_epoch);
        RV(box); RV(x); RV(x_last); RV(D); RV(D0); RV(D_last);
        RV(m_flat); RV(v_flat);
        std::uint8_t has_aux = 0; RCP_CK(has_aux);
        if (has_aux) { RV(v_update_flat); RV(m_hat_flat); RV(v_hat_flat); }
        std::uint8_t has_verlet = 0; RCP_CK(has_verlet);
        if (has_verlet)
            for (auto* nv : {&v_max, &a, &v_verlet, &a_old})
                for (auto& row : *nv) RV(row);
        RV(fixed_flags); RV(theta_anchor); RV(refresh); RV(accumulated_refresh);
        RV(particle_origin); RV(swap_D0_layout); RV(swap_log_buf); RV(spike_log_buf);
        std::uint8_t has_pp = 0; RCP_CK(has_pp);
        if (has_pp) { RV(pp_F_sq_sum); RV(pp_F_sign_flips); RV(pp_pos_min_run);
                      RV(pp_pos_max_run); RV(pp_contact_sum); RV(pp_F_sign_prev);
                      RV(pp_signflip_resets); }
        std::uint8_t has_b6a = 0; RCP_CK(has_b6a);
        if (has_b6a) RV(b6a_x_prev_step);
        RPFX(phi_history); RPFX(F_history); RPFX(U_history);
        RPFX(max_overlap_history); RPFX(mu_hist); RPFX(alpha_hist);
        RPFX(mu_flag_hist);
        RV(pair_offsets); RV(max_neighbors_per_particle);
        std::uint64_t pcap = 0; RCP_CK(pcap);
        if (ok) {
            pairs_data.reset(pcap, pair_offsets.data(), N);
            for (std::size_t i = 0; i < N; ++i) {
                std::uint32_t cnt = 0; R(&cnt, 4);
                pairs_data[pair_offsets[i]] = cnt;
                if (cnt) R(&pairs_data[pair_offsets[i] + 1],
                           static_cast<std::size_t>(cnt) * 4);
            }
        }
        std::uint8_t has_pd0 = 0; RCP_CK(has_pd0);
        if (has_pd0) RV(packed_d0_mirror);
        if (ell_on) { RV(ell_rot); RV(ell_m); RV(ell_v); RV(ell_q); RV(ell_deq); ell_sync();
                      RV(ell_cache.ids); RV(ell_cache.lam); RV(ell_cache.count); }
        if (ell_on && ell_shape.kind != 0) { RV(ell_bfac); RV(ell_rodu); RV(rod_cache.ids); RV(rod_cache.count); RV(rod_cache.tag); RV(rod_cache.nvec); }
        ok = ok && swap_engine.load_state(fp);
        char tail[4]; R(tail, 4);
        if (!ok || std::memcmp(tail, "RCPE", 4) != 0)
            throw std::runtime_error("[ckpt] truncated/corrupt checkpoint: "
                                     + config.checkpoint_path);
        std::fclose(fp);
#undef RCP_CK
        // post-restore: the epoch AND the D-stat cache are restored EXACTLY -
        // the cache is path-dependent (stale hits are part of the unbroken
        // trajectory), so bitwise resume must reproduce its staleness [B3].
        g_rcp_perm_epoch = saved_epoch;
        ++g_pairs_version;   // force the force-loop's pairs_flat mirror to resync
        start_step = static_cast<std::size_t>(st64);
        ckpt_resumed_from = start_step;
        ckpt_last_step_written = start_step;
        std::cerr << "[ckpt] resumed from step " << start_step
                  << " (" << config.checkpoint_path << ")"
                  << " swap_log=" << swap_log_buf.size() / 8
                  << " next_sweep=" << swap_next_step
                  << " interval=" << swap_cur_interval
                  << " hist=" << hist_filled << "\n";
        if (start_step + 1 > max_steps)
            throw std::runtime_error("[ckpt] checkpoint at/after max_steps; nothing to do");
    };
    if (ckpt_on) ckpt_load();

    // [101] terminal-reopen state (NOT checkpointed: a resume gets a fresh
    // reopen budget; documented semantics, keeps the v3 blob layout intact).
    long geo_reopens_done = 0;
    long geo_trades_at_reopen = -1;
    bool geo_force_step = false;
    // H12c [162]: mode-3 bneck (variant B) engages at the CONVERGED terminal via the reopen
    // hook (geo_reopens_done>0), using the shrink-to-just-touching pass. No rung-1 / skip-decay
    // machinery (variant A [160] reverted: it under-densified, bneck is not a densifier [161]).
    // The reopen "relax" is a small mu-KICK [163]: turn mu up to a small value -> slight
    // overlap (~1e-4) -> cages shift a bit -> different just-touching pairs -> a few new swaps
    // open; ADAM continues with its existing history/lr. mu_rung1 captures the rung-1 mu so the
    // kick is a small fraction of it (calibrated to keep overlaps tiny, no rearrangement).
    double mu_rung1 = 0.0;
    // [207] local quench after an injection: while step <= inj_freeze_until, every particle NOT in
    // inj_free_mask has its position and Adam moments restored after the update (only the exchange
    // participants move; kappa and the controller are untouched). Off unless swap_inject_freeze > 0.
    std::vector<std::uint8_t> inj_free_mask;
    std::size_t inj_freeze_until = 0;
    std::vector<double> frz_x, frz_m, frz_v;
    double frz_kappa = 0.0, frz_mk = 0.0, frz_vk = 0.0;
    // trial bookkeeping for accept/revert: the injected pairs, the participants' pre-exchange x/m/v, U0
    std::vector<std::pair<std::uint32_t, std::uint32_t>> inj_pairs_saved;
    std::vector<std::uint32_t> inj_participants;
    std::vector<double> inj_pre_x, inj_pre_m, inj_pre_v, inj_pre_D0;
    std::vector<std::uint64_t> inj_pre_sc;
    std::uint64_t inj_pre_geo = 0, inj_pre_cyc = 0, inj_pre_tot = 0;
    double inj_U0 = 0.0;
    // identity-keyed copies: the spatial reorder permutes every per-particle array, so index-based
    // bookkeeping does not survive a relayout; particle_origin[k] is the stable identity.
    std::vector<std::uint8_t> inj_free_origin;                         // by origin id
    std::vector<std::uint32_t> inj_part_origin;                        // participants' origin ids (participant order)
    std::vector<std::pair<std::uint32_t, std::uint32_t>> inj_pairs_origin;
    std::vector<std::uint32_t> inj_idx_of_origin;                      // scratch inverse map
    auto inj_remap_mask = [&]() {                                       // current index mask from origins
        for (std::size_t k = 0; k < N; ++k) inj_free_mask[k] = inj_free_origin[particle_origin[k]];
    };
    long inj_trials = 0, inj_reverts = 0;
    double fh_kappa_ref = kappa;   // [226] fix_height: kappa at the last box-height update

    for (std::size_t step = start_step + 1; step <= max_steps; ++step) {
        if ((step == 1 || step == 50 || step == 500) &&
                std::getenv("RCP_DEBUG_LAYOUT")) {
            std::cerr << "[rss] step" << step << ": " << rcp_rss_gb()
                      << " GB\n";
        }
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
        if (!mu_pinned && mu_decay_active && step >= mu_decay_next_tick_step && step > 250) {   // [219] suspended under mu_pin
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
        if (!mu_pinned && osc_enabled && osc_phase == OscPhase::INACTIVE &&   // [219] suspended under mu_pin
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
        if (!mu_pinned && osc_enabled && osc_phase == OscPhase::ACTIVE) {   // [219] suspended under mu_pin
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
                    // The force of the PREVIOUS step's pass (the same pass max_min_dist comes from). The slot of
                    // the current step, F_history[step - 1], is written only after this test and was always 0
                    // here, so the force condition never acted (LOG [037], BUILD.md section 10).
                    const double F_now = (step >= 2)
                        ? F_history[step - 2] : 0.0;
                    // F_now is (mean net force of the particles that carry one above 1e-8) / (mean contact force): 0/0,
                    // not a number, when no particle carries such a force, as when no contact is left. The force
                    // condition then counts as met; without this a run that had shrunk out of contact idled to the
                    // step cap, since a comparison with a NaN is false (LOG [038], [039]).
                    const bool no_contact_left = std::isnan(F_now);
                    const bool conditions_met =
                        ((F_now < stop_f_thresh) || no_contact_left) &&
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
                const double F_now = (step >= 2) ? F_history[step - 2] : 0.0;   // previous step's force: see the peak gate above
                const bool no_contact_left = std::isnan(F_now);   // 0/0: no force-carrying particle; counts as met (see the peak gate)
                const bool conditions_met =
                    ((F_now < stop_f_thresh) || no_contact_left) && (max_min_dist < osc_ov_thresh);
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
        if (!mu_pinned && osc_enabled && osc_phase == OscPhase::DONE) {   // [219] under mu_pin the run continues to max_steps
            // [101] terminal reopen (opt-in, geo only): force a geo pass and
            // re-enter FINAL_DESCENT with a fresh convergence budget. No
            // schedule constant changes; mu stays where the descent left it.
            // Stops when the budget is spent or the last reopen traded nothing.
            // H12 [148]: in terminal bottleneck mode (mode 1) the reopen pass is
            // geo_bottleneck_pass, so the sterile-stop must watch the bottleneck trade
            // counter, not geo's. Mode 2 (H12b [152]) fires bneck during the rung, not at
            // reopen, so its reopen loop watches geo trades normally.
            const long reopen_trades = (config.swap_geo_bottleneck == 1 ||
                                        config.swap_geo_bottleneck == 3)
                ? static_cast<long>(swap_engine.geo_bottleneck_trades_)
                : static_cast<long>(swap_engine.geo_trades_);
            const bool reopen_armed =
                config.swap_geo != 0 && config.swap > 0.0 &&
                config.swap_geo_reopen > 0 &&
                geo_reopens_done < config.swap_geo_reopen &&
                !(geo_reopens_done > 0 && reopen_trades == geo_trades_at_reopen);
            if (reopen_armed) {
                ++geo_reopens_done;
                geo_trades_at_reopen = reopen_trades;
                geo_force_step = true;
                osc_phase = OscPhase::FINAL_DESCENT;
                osc_consec_conv = 0;
                osc_final_start_step = step;
                mu_decay_next_tick_step = step + mu_decay_cooldown;
                // H12c [163] VARIANT B relax: at each reopen, turn mu up to a SMALL value
                // (a small fraction of the rung-1 mu) so a slight overlap (~1e-4) reopens the
                // cages a bit -> different just-touching pairs -> a few new bneck swaps. ADAM
                // keeps its history/lr. Kept tiny to avoid position rearrangement (user).
                if (config.swap_geo_bottleneck == 3 && mu_rung1 > 0.0) {
                    mu = mu_rung1 * rcp_bneck_kick_frac();
                    mu_decay_active = false;    // hold the small kick; don't immediately decay it
                }
            } else {
                steps = step;
                break;
            }
        }

        // [Cycle 19 cleanup]: per-step SEARCH and old cosine modulation
        // blocks removed (both superseded by OSC).

        // Cycle 10: periodic spatial reorder. Skips when REORDER_INTERVAL=0.
        if (REORDER_INTERVAL > 0 && step > 1 && (step % REORDER_INTERVAL == 0)) {
            reorder_particles_now();
        } else if (swap_need_relayout) {
            // SWAP capacity guard: a swapped particle outgrew the diameter its pair
            // row was provisioned for; recompute the layout (full rebuild) now.
            reorder_particles_now();
            ++swap_relayouts;
        }
        swap_need_relayout = false;


        // [219] DISTRIBUTION MORPH update - placed before D = D0*kappa so the new D0 takes effect this step.
        if (morph_on && static_cast<long>(step) >= config.morph_start && morph_lambda < 1.0 &&
            (static_cast<long>(step) - config.morph_start) % std::max<long>(1, config.morph_interval) == 0) {
            const std::size_t nD = D0.size();
            if (!morph_armed) {
                if (config.morph_target.size() != nD)
                    throw std::runtime_error("morph_target must have exactly N entries");
                morph_q_start.assign(D0.begin(), D0.end());
                std::sort(morph_q_start.begin(), morph_q_start.end());
                morph_q_target = config.morph_target;
                std::sort(morph_q_target.begin(), morph_q_target.end());
                double v0 = 0.0, v1 = 0.0;
                for (std::size_t i = 0; i < nD; ++i) {
                    v0 += morph_q_start[i] * morph_q_start[i] * morph_q_start[i];
                    v1 += morph_q_target[i] * morph_q_target[i] * morph_q_target[i];
                }
                if (!(v1 > 0.0) || morph_q_target.front() <= 0.0)
                    throw std::runtime_error("morph_target entries must be positive");
                const double sv = std::cbrt(v0 / v1);
                for (double& q : morph_q_target) q *= sv;          // volume-matched target table
                {   // a target equal to the start (to 1e-12 relative) is snapped to it exactly: a null morph is a no-op
                    double dmax = 0.0;
                    for (std::size_t i = 0; i < nD; ++i)
                        dmax = std::max(dmax, std::fabs(morph_q_target[i] / morph_q_start[i] - 1.0));
                    if (dmax < 1e-12) morph_q_target = morph_q_start;
                }
                morph_vol0 = v0;
                // band reference = max overlap at arming, floored in absolute units so a low-mu start does not
                // turn the band into an unreachably tight tolerance
                morph_ref_maxov = std::max(max_min_dist, std::max(config.morph_band_floor, 1e-7));
                morph_order.resize(nD);
                if (swap_D0_layout.size() != nD) swap_D0_layout = D0;   // capacity-guard baseline
                morph_kappa0 = kappa;
                morph_armed = true;
            }
            const double maxov_now = max_min_dist;
            // hold on either lag signal: overlaps above the band, or the global scale deflated below its arming value
            const bool held = maxov_now > config.morph_band * morph_ref_maxov ||
                              (config.morph_kappa_drop > 0.0 && kappa < morph_kappa0 * (1.0 - config.morph_kappa_drop));
            // the baseline follows kappa slowly: a quasi-static decline (the target's own phi may be lower) passes,
            // a burst (kappa falling faster than the baseline can follow) trips the guard
            morph_kappa0 += (kappa - morph_kappa0) / static_cast<double>(std::max<long>(1, config.morph_kappa_tau));
            if (held) {
                // a trip means the last growth was too fast: shrink, and lower the CEILING to the last rate that ran
                // hold-free for a full growth window (AIMD with memory -> converges to the largest sustainable rate)
                const double rmin = std::max(config.morph_rate_min, 1e-9);
                morph_rate = std::max(morph_rate * config.morph_shrink, rmin);
                morph_rate_ceiling = std::max(rmin, (morph_last_safe > 0.0) ? std::min(morph_rate_ceiling, morph_last_safe)
                                                                            : std::min(morph_rate_ceiling, morph_rate));
                morph_ok_streak = 0;
                ++morph_holds;
                if (++morph_hold_streak > std::max<long>(1, config.morph_max_hold_streak)) {
                    morph_stalled = true;          // kappa-lag stall: overlaps do not recover; stop and report
                    steps = step;
                    break;
                }
            } else {
                morph_hold_streak = 0;
                morph_lambda = std::min(1.0, morph_lambda + morph_rate);
                // grow only with a margin below the band, never above the ceiling
                if (maxov_now < config.morph_margin * config.morph_band * morph_ref_maxov) {
                    if (++morph_ok_streak >= std::max<long>(1, config.morph_grow_after)) {
                        morph_last_safe = morph_rate;                       // this rate ran a full window hold-free
                        morph_rate = std::min(morph_rate_ceiling, morph_rate * config.morph_grow);
                        morph_ok_streak = 0;
                    }
                } else {
                    morph_ok_streak = 0;
                }
                // rank-keyed apply: the particle of ascending rank r receives interpolated quantile r. The order is
                // rebuilt only when the permutation epoch moved (spatial reorder or diameter trades); the morph
                // itself preserves ranks, so at interval 1 this costs O(N), not a sort per step.
                if (morph_order_epoch != g_rcp_perm_epoch || morph_order.size() != nD) {
                    morph_order.resize(nD);
                    for (std::size_t i = 0; i < nD; ++i) morph_order[i] = static_cast<std::uint32_t>(i);
                    std::stable_sort(morph_order.begin(), morph_order.end(),
                                     [&](std::uint32_t a, std::uint32_t b) { return D0[a] < D0[b]; });
                    morph_order_epoch = g_rcp_perm_epoch;
                }
                double v = 0.0;
                bool changed = false;
                for (std::size_t r = 0; r < nD; ++r) {
                    // a + lambda*(b-a): exact identity when b == a (null morph), exact target at lambda == 1
                    const double q = (morph_lambda >= 1.0) ? morph_q_target[r]
                        : morph_q_start[r] + morph_lambda * (morph_q_target[r] - morph_q_start[r]);
                    if (q != D0[morph_order[r]]) { D0[morph_order[r]] = q; changed = true; }
                    v += q * q * q;
                }
                const double sv = std::cbrt(morph_vol0 / v);        // volume-preserving renormalisation
                if (sv != 1.0) { for (std::size_t i = 0; i < nD; ++i) D0[i] *= sv; changed = true; }
                if (changed) {   // a no-op update (null morph) must leave every derived state untouched
                    d0_cache = recompute_d0_cache();
                    D0_min_cached = d0_cache.first;
                    S_D0_cached = d0_cache.second;
                    D0_max_cached = recompute_d0_max();
                    ++g_rcp_perm_epoch;
                    morph_order_epoch = g_rcp_perm_epoch;   // ranks unchanged by the morph itself: keep the cached order
                    for (std::size_t i = 0; i < nD; ++i)
                        if (D0[i] > swap_D0_layout[i] * config.swap_relayout_ratio) { swap_need_relayout = true; break; }
                }
                ++morph_updates;
            }
            double vol_now = 0.0;
            for (std::size_t i = 0; i < nD; ++i) vol_now += D0[i] * D0[i] * D0[i];
            morph_log_buf.push_back(static_cast<double>(step));
            morph_log_buf.push_back(morph_lambda);
            morph_log_buf.push_back(morph_rate);
            morph_log_buf.push_back(maxov_now);
            morph_log_buf.push_back(held ? 1.0 : 0.0);
            morph_log_buf.push_back(vol_now);                  // sum D0^3 after the update (volume-preservation check)
        }
        // [219] settled read-out: once the target is reached, end the run when max overlap is back in band and
        // phi is flat over the settle window (the from-scratch terminal itself moves ~1e-3 per 500 steps in its
        // oscillation phase, so the tolerance is a user choice; 0 = run to max_steps).
        if (morph_on && morph_armed && morph_lambda >= 1.0 && config.morph_settle_flat > 0.0 &&
            static_cast<long>(step) > config.morph_settle_window + 2 &&
            max_min_dist <= config.morph_band * morph_ref_maxov) {
            const std::size_t w = static_cast<std::size_t>(std::max<long>(2, config.morph_settle_window));
            auto lo = phi_history.begin() + (step - 1 - w), hi = phi_history.begin() + (step - 1);
            auto mm = std::minmax_element(lo, hi);
            if (*mm.second - *mm.first < config.morph_settle_flat) {
                steps = step;
                break;
            }
        }

        g_t_bk_alpha.begin();
        // phi_history holds phi / phi_modifier, so the comparison is made with phi / phi_modifier. The released code
        // compared the raw phi: in a curved container (phi_modifier != 1) "phi fell by more than 5e-4" was then always
        // true and alpha was cut by 1.5 every 500 steps for the whole run (LOG [035], [040]; repaired for every shape).
        const double phi_cmp = phi / phi_modifier;
        if ((step - LastAlphaUpdate > 500) &&
            (method == "ADAM") &&
            ((phi_cmp - phi_history[std::max<std::size_t>(step - 250, 1)]) < 5e-4) &&
            (step > 2500) && step > mu_change + 500) {

            double trend1 = mean(F_history, step - 450, step - 250) / mean(F_history, step - 75, step - 1);
            double trend2 = (phi_cmp - phi_history[std::max<std::size_t>(step - 250, 1)]);

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
            // Cycle 21k: parallel — each D[i] is one rounding of
            // D0[i]*kappa regardless of thread, so values are
            // byte-identical to the serial scan.
            const std::size_t nD = D0.size();
            if (!has_fixed) {
                #pragma omp parallel for schedule(static)
                for (std::size_t i = 0; i < nD; ++i) {
                    D[i] = D0[i] * kappa;
                }
            } else {
                // Cascade: frozen particles keep their diameter (D0); only the
                // active class scales with kappa. Branch hoisted out of the loop.
                #pragma omp parallel for schedule(static)
                for (std::size_t i = 0; i < nD; ++i) {
                    D[i] = fixed_flags[i] ? D0[i] : D0[i] * kappa;
                }
            }
        }
        g_t_bk_dupdate.end();
        // Cycle 4a: O(1) replacements for the two O(N) D-derived scans below.
        // Dmin_lastest = min(D), but D[i] = D0[i] * kappa for all i, so
        //   min(D) = D0_min_cached * kappa.
        // When runtime_fix_diameter is true, D == D0 and kappa is held at 1
        // (see target_phi branch above), so D0_min_cached IS Dmin.
        double Dmin_lastest = has_fixed
            ? *std::min_element(D.begin(), D.end())          // cascade: mixed scaling
            : D0_min_cached * (runtime_fix_diameter ? 1.0 : kappa);

        coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
        // current_volume = sum_i coeff * (D[i]/2)^Ndim
        //                = coeff * S_D0_cached * kappa^Ndim   (or kappa^Ndim = 1
        //                  when runtime_fix_diameter and kappa is locked at 1).
        double kappa_to_Ndim = 1.0;
        if (!runtime_fix_diameter) {
            for (std::size_t k = 0; k < Ndim; ++k) kappa_to_Ndim *= kappa;
        }
        if (has_fixed) {
            // cascade: frozen particles do not scale with kappa, so the O(1)
            // S_D0*kappa^Ndim shortcut is invalid; recompute the volume directly.
            current_volume = 0.0;
            for (std::size_t i = 0; i < N; ++i)
                current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
        } else {
            current_volume = coeff * S_D0_cached * kappa_to_Ndim;
        }

        if (fix_height && !runtime_fix_diameter && kappa != fh_kappa_ref) {
            // [226] the height follows the particle scale: H = H_start x kappa. Scale by the CHANGE of kappa since
            // the last update (scaling by kappa itself every step compounded to kappa^steps and overflowed).
            const double r = kappa / fh_kappa_ref;
            box[Ndim - 1] = box[Ndim - 1] * r;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i * Ndim + (Ndim - 1)] = x[i * Ndim + (Ndim - 1)] * r;
            }
            fh_kappa_ref = kappa;
        }

        box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());

        phi = current_volume / box_volume;

        if (options.has_target_phi && !target_phi_locked && !runtime_fix_diameter) {
            const double user_phi = phi / phi_modifier;
            if (user_phi >= options.target_phi) {
                auto target_scaled = scale_diameters_nd(
                    D,
                    ell_on ? options.target_phi * phi_modifier / ell_vf : options.target_phi * phi_modifier,
                    box,
                    Ndim,
                    fix_height);
                D = std::move(target_scaled.first);
                if (fix_height) {   // [226] the height follows the diameters' rescale at the lock
                    box[Ndim - 1] *= target_scaled.second;
                    for (std::size_t i = 0; i < N; ++i) x[i * Ndim + (Ndim - 1)] *= target_scaled.second;
                    box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());
                }
                D0 = D;
                kappa = 1.0;
                fh_kappa_ref = 1.0;
                mu = 0.0;
                runtime_fix_diameter = true;
                target_phi_locked = true;

                // D0 has been rescaled — refresh the cached quantities.
                d0_cache = recompute_d0_cache();
                D0_min_cached = d0_cache.first;
                S_D0_cached   = d0_cache.second;
                D0_max_cached = recompute_d0_max();   // Cycle 21k

                Dmin_lastest = D0_min_cached;  // kappa = 1 here
                coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
                current_volume = coeff * S_D0_cached;  // kappa^Ndim = 1
                phi = current_volume / box_volume;
                // Cycle 13: D was rescaled, so per-particle MaxNb caps shift.
                {
                    std::uint64_t new_total = compute_pair_layout(
                        max_neighbors_per_particle, pair_offsets);
                    pairs_data.reset(new_total, pair_offsets.data(), N);
                }
                std::fill(refresh.begin(), refresh.end(), 1u);
                {
            std::uint32_t call_max = 0, call_min = std::numeric_limits<std::uint32_t>::max();
            get_pairs_dispatch(N, Ndim, x.data(), D, box, walls, refresh,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(),
                               pairs_data.data(),
                               &call_max, &call_min,
                               packed_d0_ptr(), D0.data());
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
                x_last = x;
                D_last = D;
                Dmin = Dmin_lastest;
                shell_refill_anchor();   // Cycle 21k
            }
        }

        phi_history[step - 1] = phi / phi_modifier;

        if (phi > phi_max.first) {
            phi_max = std::make_pair(phi, step);
        }

        if ((step > 500 && mu_flag == 1)) {
            // Cycle 23: phi-plateau early-out for the rung 1->0 transition
            // DISABLED (see trigger below). XX/delta_XX computed the phi range
            // over the last 250 steps; no longer used.
            // std::vector<double> XX(phi_history.begin() + step - 250, phi_history.begin() + step - 1);
            // double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

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
                // Cycle 21M: absolute -2000 cut, all presets (Ken's spec).
                // Safe only WITH the default-on rung-1 mu boost, which
                // front-loads compression so the shortened schedule still
                // reaches the over-compression state the oscillation
                // phase needs (3-seed dose-response at 100K: deeper cuts
                // lose phi even with the boost — the knee is sharp).
                // The oscillation period follows automatically
                // (osc_period = 0.5 * rung-1 end step). Floor keeps
                // tiny-N immediate-preset runs from degenerating.
                double dl = static_cast<double>(Y) * 4000.0 + s_term
                            - 2000.0;
                if (dl < 1500.0) dl = 1500.0;
                return static_cast<std::size_t>(dl);
            }();
            // Cycle 23: rung 1->0 transition gated on the rung-1 DEADLINE ONLY.
            // The two early-outs below (phi-plateau delta_XX<5e-6, and no-new-
            // phi-max-in-3500-steps stall) fired prematurely at large N for
            // narrow distributions (Weibull k>=1, and analogous power-law cases),
            // ending rung-1 early AND shortening the downstream oscillation
            // period (osc = 0.5 * rung-1 end step) -> under-compaction, phi
            // ~0.003-0.004 low and non-monotonic vs Dmax/L. Safety fix: always
            // run the full rung-1 budget. (Costs wall time on fast-converging
            // cases; revisit later.) Original OR-condition preserved here:
            // if (delta_XX < 5e-6 || step - phi_max.second > 3500 || (step > static_cast<int>(rung1_deadline))) {
            const std::size_t mu_hold_shift =
                (config.mu_hold_shift_deadline && config.mu_hold_steps > 0 &&
                 config.mu_hold_start >= 0 &&
                 static_cast<std::size_t>(config.mu_hold_start) <= rung1_deadline)
                ? static_cast<std::size_t>(config.mu_hold_steps) : 0;
            if (!mu_pinned && step > static_cast<int>(rung1_deadline + mu_hold_shift)) {   // [219] suspended under mu_pin
                // Cycle 19 Step 2: arm adaptive feedback mu decay.
                // mu / 10 (= mu_1 / 10): empirically the winning ratio.
                // Keeps mu_1 → mu_0 a decade decay, so t_relax measured
                // here is consistent with the decade decay used in the
                // oscillation reference exp(-t/t_relax).
                // [160->162] VARIANT A (skip-decay bneck-convergence) REVERTED: it
                // under-densified (bneck is a permutation, not a densifier [161]). The full
                // normal convergence runs for ALL modes; mode-3 bneck now engages at the
                // CONVERGED terminal via the reopen hook below (variant B, user [latest]).
                mu_decay_target = mu / 10.0;
                mu_decay_active = true;
                mu_decay_next_tick_step = step + mu_decay_cooldown;
                alpha /= 2.0;
                mu_flag = 0;
                mu_change = step;
                if (config.swap_geo_bottleneck == 3 && mu_rung1 <= 0.0) mu_rung1 = mu;   // H12c [163]: capture rung-1 mu for the reopen kick (mode 3 only)
                // Cascade (additive, has_fixed-guarded): release the frozen-
                // diameter scaffold at the FIRST rung-1->0 trigger so every class
                // is now free to densify / shrink-to-contact together (the McGeary
                // final settle). Re-anchor D0 for the frozen particles so that
                // D = D0*kappa stays exactly continuous across the release (no size
                // jump), then clear the mask. has_fixed stays true, so Dmin /
                // current_volume keep their O(N) recompute (no stale cache from the
                // changed D0), and the all-zero mask makes the D-update and the
                // dkappa/pressure terms reduce identically to the standard path.
                if (has_fixed && !cascade_released) {
                    for (std::size_t ci = 0; ci < fixed_flags.size(); ++ci) {
                        if (fixed_flags[ci]) D0[ci] = D[ci] / kappa;
                        fixed_flags[ci] = 0;
                    }
                    cascade_released = true;
                }
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
        if (!mu_pinned && !osc_enabled && step > 500 && mu_flag == 0 && !mu_decay_active &&   // [219]
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
                               &call_max, &call_min,
                               packed_d0_ptr(), D0.data());
            if (call_max > agg_max_neighbors) agg_max_neighbors = call_max;
            if (call_min < agg_min_neighbors) agg_min_neighbors = call_min;
        }
            g_t_pairs.end();
            g_t_bookkeep.begin();
            x_last = x;
            D_last = D;
            Dmin = Dmin_lastest;
            shell_refill_anchor();   // Cycle 21k
        } else {
            std::fill(refresh.begin(), refresh.end(), 0u);
            g_t_bookkeep.end();
            g_t_refresh.begin();
            // Cycle 21k adaptive shell controller. Signal: EMA of the
            // per-step kappa growth rate. Worst-case per-step budget burn
            // is on the largest particle (dD/2 = Dmax*gr/2); choose s so a
            // freshly-anchored budget survives ~RCP_SHELL_K steps of pure
            // growth, capped at RCP_SHELL_SMAX. Steady state (gr ~ 0)
            // decays s back to 1 and reproduces fixed-shell behavior.
            {
                // SIGNED growth rate: kappa oscillates mid-run (dips and
                // recoveries), and a positive-only signal rectifies that
                // oscillation into phantom growth, keeping the shell wide
                // for the whole run (measured: +41% force-loop time). The
                // signed EMA lets dips cancel recoveries, so only genuine
                // net compression (the early phase) widens the shell.
                const double gr_step = (shell_kappa_prev > 0.0)
                    ? (kappa / shell_kappa_prev - 1.0) : 0.0;
                shell_kappa_prev = kappa;
                // The EMA may go negative during dip phases — that is the
                // point (oscillation nets out); the burn>0 guard below is
                // the only rectifier.
                shell_gr_ema = 0.9 * shell_gr_ema + 0.1 * gr_step;
                double s_dyn = 1.0;
                const double smax = rcp_shell_smax();
                if (smax > 1.0 && !runtime_fix_diameter) {
                    const double theta_base =
                        Dmin * 0.25 * trigger_scale_run;
                    const double burn =
                        0.5 * D0_max_cached * kappa * shell_gr_ema;
                    if (burn > 0.0 && theta_base > 0.0) {
                        const double s_need =
                            rcp_shell_k() * burn / theta_base;
                        s_dyn = std::min(std::max(1.0, s_need), smax);
                    }
                }
                g_rcp_shell_scale_dyn = s_dyn;
            }
            // Cycle 20 L2: trigger scale resolved once per run (see
            // trigger_scale_run above for the safety derivation).
            const double dmin_quarter = Dmin * 0.25 * trigger_scale_run;
            // Budget granted to particles re-anchored THIS step (their
            // lists rebuild under the current shell scale).
            const double theta_new_anchor =
                dmin_quarter * g_rcp_shell_scale_dyn;
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
            // Cycle 20 growth-aware trigger: diameter growth since the
            // anchor consumes the same pair-shell budget displacement does
            // — the pair cutoff (Di+Dj)/2 grows while the stale list's
            // membership was decided with the old cutoff. The sentinel
            // caught real missed pairs from this on BOTH backends in the 2D
            // early phase.
            //
            // Cycle 21j: kappa OSCILLATES (the ADAM kappa update overshoots
            // during compression), so growth since the anchor is NOT
            // monotone. The original max(0, D - D_anchor) clamp let
            // recovery growth after a kappa dip go uncounted: a pair list
            // built near the dip expands its cutoff budget-free on the way
            // back up (sentinel-caught: negative dD/2 on the partner of a
            // missed pair). Fix: D_last[k] now tracks the RUNNING MINIMUM
            // of D[k] since the anchor, so growth = D - minD >= 0 bounds
            // the cutoff expansion relative to ANY possible list-build
            // moment in the anchor window. For monotone growth this is
            // identical to the original anchor semantics.
            #pragma omp parallel for reduction(||:update_flag) schedule(static)
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t base = k * Ndim;
                double d2 = 0.0;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double dd = x[base + d] - x_last[base + d];
                    d2 += dd * dd;
                }
                if (D[k] < D_last[k]) D_last[k] = D[k];   // running min
                const double growth_half = (D[k] - D_last[k]) * 0.5;
                // Cycle 21k: per-particle budget, min-clamped to the
                // CURRENT scale's budget. The clamp makes the bound hold
                // for any s(t) trajectory: a stale particle's budget never
                // exceeds theta at any partner's later rebuild time, so
                // theta_j + 2*theta_k <= 3*theta(s_build) <= t_shell(s_build).
                if (theta_anchor[k] > theta_new_anchor)
                    theta_anchor[k] = theta_new_anchor;
                const double budget = theta_anchor[k] - growth_half;
                if (budget <= 0.0 || d2 > budget * budget) {
                    refresh[k] = 1;
                    for (std::size_t d = 0; d < Ndim; ++d)
                        x_last[base + d] = x[base + d];
                    D_last[k] = D[k];
                    theta_anchor[k] = theta_new_anchor;   // re-anchor
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
                #pragma omp parallel for schedule(static)
                for (std::size_t k = 0; k < N; ++k) {
                    D_last[k] = D[k];
                    theta_anchor[k] = theta_new_anchor;   // Cycle 21k
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
                               &call_max, &call_min,
                               packed_d0_ptr(), D0.data());
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
        // Cycle 21M: rung-1 mu boost (no-op at B=1; mu itself untouched).
        const bool mu_hold_now = config.mu_hold_steps > 0 && config.mu_hold_start >= 0 &&
            static_cast<long>(step) >= config.mu_hold_start &&
            static_cast<long>(step) < config.mu_hold_start + config.mu_hold_steps;
        if (mu_hold_now && config.swap_mu_ctrl && config.mu_ctrl_interval > 0 &&
            (static_cast<long>(step) - config.mu_hold_start) % config.mu_ctrl_interval == 0) {
            std::size_t nz = 0;
            for (std::size_t q = 0; q < N; ++q) nz += (z[q] >= 1);
            const double v = static_cast<double>(nz) / static_cast<double>(N);
            const double err = config.mu_ctrl_v_target - v;
            if (std::fabs(err) > config.mu_ctrl_deadband)
                mu_ctrl_mult *= std::exp(config.mu_ctrl_eta * err);
            mu_ctrl_mult = std::min(config.mu_ctrl_max_mult, std::max(config.mu_ctrl_min_mult, mu_ctrl_mult));
        }
        const double mu_hold_mult_now = config.swap_mu_ctrl ? mu_ctrl_mult : config.mu_hold_mult;
        const double mu_eff = mu_pinned ? config.mu_pin :             // [219] absolute mu pin
            mu_hold_now ? mu_hold_mult_now * mu1_ref :
            (mu_flag == 1 && rcp_mu1_boost() > 1.0)
            ? mu * (1.0 + (rcp_mu1_boost() - 1.0) *
                          std::exp(-static_cast<double>(step) /
                                   rcp_mu1_boost_tau()))
            : mu;
        if (mu_hold_now) ++mu_hold_applied;
        // Cascade: per-particle kappa-radius (fixed -> 0, excluded from dkappa) and
        // the active-only diameter sum for the pressure term. Built only when
        // has_fixed; otherwise pass nullptr/-1 -> get_forces_nd_3 byte-identical.
        const double* rk_ptr = nullptr; double sumD_active = -1.0;
        std::vector<double> rk_flat_cascade;
        if (has_fixed) {
            rk_flat_cascade.resize(N); sumD_active = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                if (fixed_flags[i]) { rk_flat_cascade[i] = 0.0; }
                else { rk_flat_cascade[i] = D[i] * 0.5; sumD_active += D[i]; }
            }
            rk_ptr = rk_flat_cascade.data();
        }
        // SWAP: squared pair separations are stashed by the force loop (only when
        // swap is on) so the sweep below needs no second pass over the pairs.
        // SWAP: on a sweep step, the force loop collects the near pairs the sweep
        // needs (one compare per pair, no second pass over the pair lists).
        const bool swap_step_now =
            swap_on && !config.swap_sweeps_off &&
            !runtime_fix_diameter && step >= swap_next_step &&
            (swap_params.stop < 0 || static_cast<long>(step) <= swap_params.stop) &&
            (swap_params.phase == 0 ||
             (mu_flag == 1 && (swap_params.phase & 1)) ||
             (mu_flag == 0 && (swap_params.phase & 2)) ||
             (mu_flag == -1 && (swap_params.phase & 4)));
        // rescue channel: collect near-pair edges (for exact dE) only every
        // rescue_interval-th sweep step, so the collection cost is paid rarely
        // [207] while a trial quench is active the regular swap passes are suspended: only the
        // injected exchange is under test, and a revert must land restored positions on the
        // pre-exchange diameters (four geo passes in a 20-step quench had re-traded ~all participants).
        const bool trial_quench_active = (step <= inj_freeze_until + 1) && !inj_free_mask.empty();   // +1: through the decision step
        const bool rescue_step_now =
            swap_step_now && !trial_quench_active &&
            (swap_params.rescue > 0.0 || swap_params.cycle3 > 0.0) &&
            config.swap_exact == 0 &&
            (swap_params.rescue_interval <= 1 ||
             step % static_cast<std::size_t>(swap_params.rescue_interval) == 0);
        const bool geo_step_now =
            config.swap_geo != 0 && config.swap > 0.0 && !trial_quench_active &&
            !quench_on &&                      // [225] a quench runs at frozen diameters: no permutation passes
            (geo_force_step ||                 // [101] reopen-forced pass, bypasses window/cadence
             (static_cast<long>(step) >= config.swap_geo_start &&
              (config.swap_geo_stop < 0 || static_cast<long>(step) < config.swap_geo_stop) &&
              ((static_cast<long>(step) - config.swap_geo_start) % std::max<long>(1, config.swap_geo_k) == 0)));
        const bool is_reopen_pass = geo_force_step;   // [172] capture BEFORE reset (used at the reopen geo_pass for nonsep==2)
        if (geo_step_now) geo_force_step = false;
        const bool sweep_step_now = swap_step_now && !trial_quench_active && !(geo_step_now && config.swap_geo == 1);
        if (quench_on && !(swap_step_now || geo_step_now)) { swap_sink.want_edges = false; swap_sink.want_gsum = true; }   // [221]
        if (swap_step_now || geo_step_now) {
            swap_sink.want_edges = (config.swap_exact != 0) || rescue_step_now || geo_step_now ||
                                   (config.swap_event_log_geom != 0 && config.swap_event_log_stride > 0);
            swap_sink.want_gsum = (config.swap_exact == 0);
            if (ell_on) swap_sink.want_edges = false;   // the ellipsoid pass reads the neighbour lists itself
            if (swap_sink.want_edges) swap_sink.margin = swap_engine.edge_margin(D, D0);
            if ((config.swap_event_log_geom != 0 && config.swap_event_log_stride > 0) || geo_step_now) {
                // ME0-iv: the geometry shell is r < 1.1 (R_i + R_k); the geo pass uses
                // the 1.5x shell -> cover r - (R_i+R_k) <= 0.5 D_max on geo steps
                double dmax = 0.0;
                for (std::size_t q = 0; q < N; ++q) if (D[q] > dmax) dmax = D[q];
                swap_sink.margin = std::max(swap_sink.margin, (geo_step_now ? 0.50 : 0.10) * dmax);
            }
        }
        if (ell_on) {
            shape_forces(pairs_data.data(), pair_offsets.data(), x.data(), N, ell_shape, D, box, walls, ell_rot.data(),
                       F, ell_torque.data(), U, max_min_dist, Lc, Fmean, mu_eff, dkappa, z, ell_energy, ell_amat,
                       ell_shape.track_tol > 0.0 ? &ell_cache : nullptr,
                       quench_on ? (swap_sink.gsum.resize(N), swap_sink.gsum.data()) : nullptr,
                       rod_cache.count.empty() ? nullptr : &rod_cache);
            if (ell_capture_armed) {   // test hook: the engine's own numbers on the state right after exchanges
                ell_capture_armed = false;
                const std::size_t nrot = ell_nrot(Ndim);
                ell_cap_step = static_cast<long>(step); ell_cap_energy = ell_energy;
                ell_cap_x = x; ell_cap_D = D; ell_cap_rot = ell_rot;
                ell_cap_q.assign(N * (Ndim - 1), 1.0);
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t k = 1; k < Ndim; ++k) ell_cap_q[i * (Ndim - 1) + k - 1] = ell_q[i * Ndim + k];
                ell_cap_F.assign(F, F + N * Ndim);
                ell_cap_T.assign(ell_torque.begin(), ell_torque.begin() + N * nrot);
                ell_cap_z.assign(N, 0.0);
                for (std::size_t i = 0; i < N; ++i) ell_cap_z[i] = static_cast<double>(z[i]);
            }
        } else
        get_forces_nd_3(pairs_data.data(), pair_offsets.data(), x.data(), N, Ndim, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu_eff, dkappa, z,
                        // Cycle 21k: packed-D0 + the kappa the last D-update used
                        // (D[j] == D0[j]*kappa bit-for-bit at this point).
                        packed_d0_on ? packed_d0_mirror.data() : nullptr, kappa,
                        rk_ptr, sumD_active, (swap_step_now || geo_step_now || quench_on) ? &swap_sink : nullptr);
        g_t_forces.end();

        // ---- diameter-SWAP sweep (additive; skipped entirely when swap == 0) ------
        // Runs right after the force evaluation: x, D and the pair lists are the
        // ones the forces were just computed with, and swap_pair_d2 holds their
        // exact squared separations. Accepted swaps exchange D0/D of two particles
        // at fixed positions; the position update below still uses the forces of
        // the pre-swap diameters (one-step lag, as with any Adam step), and the next
        // step's D-update / refresh check see the exchanged D0. Neighbour-list
        // validity: growth of a swapped particle is charged to its refresh budget
        // through the running-min D_last (growth_half), exactly like kappa growth,
        // so no forced refresh is needed for small |dD|. The kd-tree D-stat index
        // cache is invalidated (perm epoch) because an exchange leaves the D[0] and
        // D[N/2] sentinels intact while moving the max/min/median indices.
        swap_engine.current_step_ = step;
        if (swap_step_now || geo_step_now) {
            g_t_swap.begin();
            swap_touched.clear();
            SwapSweepStats sst;
            const bool geo_sphere_now = geo_step_now && !ell_on;   // the sphere engine's passes below
            if (geo_step_now && ell_on) {
                // ---- ellipsoids: one permutation pass (ellipsoid.hpp, ell_permute_pass). Two particles trade
                // places: diameter, shape and (ell_permute_orient) orientation, at fixed positions.
                EllPermuteParams ep;
                ep.window = config.swap_geo_w;
                ep.budget = std::max<long>(0, config.swap_geo_B);
                ep.max_partners = config.ell_permute_partners;
                ep.move_orientations = config.ell_permute_orient != 0;
                ep.decide_in_place = config.ell_permute_orient != 1;
                ep.overlap_guard = config.ell_permute_guard != 0;
                ep.pass_index = swap_engine.geo_passes_;
                ep.n_iter = std::max(1, config.ell_iters);
                ep.kind = ell_shape.kind;
                ell_exch.clear();
                ell_perm_candidates += ell_permute_pass(pairs_data.data(), pair_offsets.data(), x.data(), N, Ndim, box,
                                                        walls, D.data(), D0.data(), ell_q.data(), ell_deq.data(),
                                                        ell_rot.data(), ell_amat.data(), ep, ell_exch);
                ++swap_engine.geo_passes_;
                swap_engine.geo_trades_ += ell_exch.size();
                if (!ell_exch.empty()) {
                    auto& sc = swap_engine.swap_counts();
                    if (sc.size() != N) sc.assign(N, 0);
                    const std::size_t nrot = ell_nrot(Ndim);
                    std::vector<char> tch(N, 0);
                    for (const EllExchangeRecord& e : ell_exch) {
                        swap_touched.push_back(e.i); swap_touched.push_back(e.j);
                        std::swap(ell_bfac[e.i], ell_bfac[e.j]);
                        if (!ell_rodu.empty()) { std::swap(ell_rodu[2 * e.i], ell_rodu[2 * e.j]); std::swap(ell_rodu[2 * e.i + 1], ell_rodu[2 * e.j + 1]); }
                        ++sc[e.i]; ++sc[e.j]; tch[e.i] = tch[e.j] = 1;
                        if (ep.move_orientations) {
                            // the orientation's Adam state: the second moment (the scale of the torques this
                            // particle has felt) travels with it, the first moment (a direction that belonged
                            // to the old site) is dropped
                            for (std::size_t c = 0; c < nrot; ++c) {
                                std::swap(ell_v[e.i * nrot + c], ell_v[e.j * nrot + c]);
                                ell_m[e.i * nrot + c] = 0.0; ell_m[e.j * nrot + c] = 0.0;
                            }
                        }
                    }
                    // a pair with an exchanged particle starts its contact function cold on the next pass
                    for (std::size_t i = 0; i < N; ++i) {
                        if (tch[i]) { ell_cache.count[i] = 0; continue; }
                        std::uint32_t keep = 0;
                        for (std::uint32_t c = 0; c < ell_cache.count[i]; ++c) {
                            const std::uint32_t id = ell_cache.ids[i * EllPairCache::CAP + c];
                            if (tch[id]) continue;
                            ell_cache.ids[i * EllPairCache::CAP + keep] = id;
                            ell_cache.lam[i * EllPairCache::CAP + keep] = ell_cache.lam[i * EllPairCache::CAP + c];
                            ++keep;
                        }
                        ell_cache.count[i] = keep;
                    }
                    if (!rod_cache.count.empty()) rod_cache.purge(tch);
                    sst.accepted = ell_exch.size();
                    ++ell_perm_productive;
                    if (config.ell_capture > 0 && ell_perm_productive == config.ell_capture) {
                        ell_capture_armed = true;
                        ell_cap_exch.clear();
                        for (const EllExchangeRecord& e : ell_exch) {
                            ell_cap_exch.push_back(static_cast<double>(e.i)); ell_cap_exch.push_back(static_cast<double>(e.j)); ell_cap_exch.push_back(e.gain);
                        }
                    }
                }
            }
            if (sweep_step_now) sst = swap_engine.sweep(N, D, D0, swap_sink, z, swap_touched);
            // H7 HANDOFF [127]: when capk is on with a handoff step, geo runs ONLY
            // before it and cage runs ONLY after it -> no competing every-step passes
            // (fixes the [126] thrash). capk_start==0 keeps the old behavior (both).
            const bool capk_on = (config.swap_geo_capk != 0);
            const long capk_start = config.swap_geo_capk_start;
            const bool cage_on = (config.swap_geo_cage != 0);       // H8 [131]
            const long cage_start = config.swap_geo_cage_start;
            // geo yields once EITHER handoff (capk or cage) has fired: c_min-geo seeds
            // the state for the first M steps, then the gated pass mops up the residual.
            const bool geo_yielded =
                (capk_on && capk_start > 0 && static_cast<long>(step) >= capk_start) ||
                (cage_on && cage_start > 0 && static_cast<long>(step) >= cage_start);
            if (geo_sphere_now && !geo_yielded) {
                // H12 [148]/[162]: at a reopen (= the CONVERGED terminal, OSC-DONE) the forced
                // pass is bottleneck-relief; normal c_min-geo runs during initial compression.
                // mode 1 = un-shrunk bneck (infidelic [155]); mode 3 = shrink-to-just-touching
                // bneck (variant B, engages at the converged terminal 0.8806, user [latest];
                // at worst 0 swaps -> phi unchanged). Mode 2 leaves geo_pass here (superseded).
                if ((config.swap_geo_bottleneck == 1 || config.swap_geo_bottleneck == 3)
                    && geo_reopens_done > 0)
                    swap_engine.geo_bottleneck_pass(N, D, D0, swap_sink, swap_touched);
                else {
                    swap_engine.geo_reopen_active_ = is_reopen_pass;   // [172] H19: coupled rescore only at reopen when geo_nonsep==2
                    swap_engine.geo_pass(N, D, D0, swap_sink, swap_touched);
                    swap_engine.geo_reopen_active_ = false;
                }
            }
            if (geo_sphere_now && config.swap_geo_cycle != 0 && !geo_yielded)  // H2 [117]: 3-cycles
                swap_engine.geo_cycle_pass(N, D, D0, swap_sink, swap_touched);
            if (geo_sphere_now && config.swap_geo_cycle4 != 0 && !geo_yielded) // H34 [194]: 4-cycles
                swap_engine.geo_cycle4_pass(N, D, D0, swap_sink, swap_touched);
            // H12b [152]: mode 2 = bottleneck-relief ADDITIVELY during the mu_1 rung
            // (mu_flag==1), under residual overlap/compression above jamming ("the signal
            // is the overlap"). geo_start already gates to the dense post-onset regime;
            // the rung-1 ADAM dynamics provide the gentle relax between passes.
            // NOTE [155]: mode 2 evaluates at the CURRENT (overlapping) scale and does NOT
            // shrink -> starved accommodation gate -> 0 trades (infidelic). Superseded by 3.
            if (geo_sphere_now && config.swap_geo_bottleneck == 2 && mu_flag == 1 && !geo_yielded)
                swap_engine.geo_bottleneck_pass(N, D, D0, swap_sink, swap_touched);
            // H12c [162] VARIANT B: mode-3 bneck engages at the CONVERGED terminal (not
            // during the run) via the reopen hook (geo_force_step), where geo_reopens_done>0.
            // Handled in the forced-pass branch above (config.swap_geo_bottleneck 1||3 ->
            // geo_bottleneck_pass). At worst it makes 0 swaps and phi is unchanged; it can
            // only keep or improve phi, never lose density.
            if (geo_sphere_now && capk_on &&
                (capk_start <= 0 || static_cast<long>(step) >= capk_start))  // H7 cage (post-handoff)
                swap_engine.geo_capk_pass(N, D, D0, swap_sink, swap_touched);
            if (geo_sphere_now && cage_on &&
                (cage_start <= 0 || static_cast<long>(step) >= cage_start))  // H8 overlap-energy (post-handoff)
                swap_engine.geo_cage_pass(N, D, D0, swap_sink, swap_touched,
                                          x.data(), box.data(), static_cast<int>(Ndim));  // x/box used only if geo_cage_move [133]
            if (!swap_touched.empty()) {
                ++g_rcp_perm_epoch;
                for (std::uint32_t k : swap_touched) {
                    if (config.swap_refresh) theta_anchor[k] = 0.0;   // optional forced refresh
                    if (D0[k] > swap_D0_layout[k] * config.swap_relayout_ratio) swap_need_relayout = true;
                    if (swap_params.adam >= 1) {
                        for (std::size_t d = 0; d < Ndim; ++d) m_flat[k * Ndim + d] = 0.0;
                    }
                    if (swap_params.adam >= 2) {
                        for (std::size_t d = 0; d < Ndim; ++d) v_flat[k * Ndim + d] = 0.0;
                    }
                }
                swap_total_accepted += sst.accepted;
            }
            // Adaptive back-off: an idle sweep (accepted <= idle_frac * N -- relative,
            // so that a handful of stragglers at large N do not keep every step
            // sweeping) doubles the interval (up to interval*backoff); a productive
            // sweep resets it. Cost follows activity.
            const bool sweep_idle = config.swap_ctrl
                ? sst.starved
                : (static_cast<double>(sst.accepted) <=
                   config.swap_idle_frac * static_cast<double>(N));
            if (!sweep_step_now) {
                // geo-only step: leave the sweep cadence and its log untouched
            } else
            if (sweep_idle && config.swap_backoff > 1) {
                const std::size_t cap = static_cast<std::size_t>(swap_params.interval) *
                                        static_cast<std::size_t>(config.swap_backoff);
                swap_cur_interval = std::min(cap, swap_cur_interval * 2);
            } else {
                swap_cur_interval = static_cast<std::size_t>(swap_params.interval);
            }
            if (sweep_step_now) swap_next_step = step + swap_cur_interval;
            if (sweep_step_now) {
            swap_log_buf.push_back(static_cast<double>(step));
            swap_log_buf.push_back(static_cast<double>(mu_flag));
            swap_log_buf.push_back(static_cast<double>(sst.candidates));
            swap_log_buf.push_back(static_cast<double>(sst.accepted));
            swap_log_buf.push_back(sst.dE_sum);
            swap_log_buf.push_back(phi / phi_modifier);
            swap_log_buf.push_back(static_cast<double>(sst.contacts));
            // slot 8: mean contact energy (exact mode) or the controller's window
            swap_log_buf.push_back(config.swap_ctrl ? sst.w_used : sst.e_scale);
            }
            g_t_swap.end();
        }
        // ---- surgical injections (ME0-ii [log 054]) ---------------------------
        if (swap_on && !config.swap_inject.empty()) {
            static thread_local std::vector<std::uint32_t> inj_touched;
            inj_touched.clear();
            std::size_t n_inj = 0;
            static thread_local std::vector<double> inj_D0_before;            // [207] pre-injection D0 (for the revert check)
            if (config.swap_inject_freeze > 0) inj_D0_before = D0;
            for (std::size_t q = 0; q + 2 < config.swap_inject.size(); q += 3) {
                if (static_cast<std::size_t>(config.swap_inject[q]) != step) continue;
                if (swap_engine.inject(N, D, D0, inj_touched,
                        static_cast<std::uint32_t>(config.swap_inject[q + 1]),
                        static_cast<std::uint32_t>(config.swap_inject[q + 2])))
                    ++n_inj;
            }
            if (!inj_touched.empty()) {
                ++g_rcp_perm_epoch;
                // [206] injection-only Adam reset (a neighbourhood exchange invalidates the moments of
                // its participants; ordinary swaps keep swap_adam). -1 = follow swap_adam (legacy).
                const int inj_adam = config.swap_inject_adam >= 0 ? config.swap_inject_adam : swap_params.adam;
                for (std::uint32_t k : inj_touched) {
                    if (config.swap_refresh) theta_anchor[k] = 0.0;
                    if (D0[k] > swap_D0_layout[k] * config.swap_relayout_ratio)
                        swap_need_relayout = true;
                    if (inj_adam >= 1)
                        for (std::size_t d = 0; d < Ndim; ++d) m_flat[k * Ndim + d] = 0.0;
                    if (inj_adam >= 2)
                        for (std::size_t d = 0; d < Ndim; ++d) v_flat[k * Ndim + d] = 0.0;
                }
                swap_total_accepted += n_inj;
                if (config.swap_inject_freeze > 0 && n_inj > 0) {          // [207] arm the local quench
                    inj_free_mask.assign(N, 0);
                    for (std::uint32_t k : inj_touched) inj_free_mask[k] = 1;
                    // user [207]: the participants alone cannot force-balance against a rigid wall -
                    // free their neighbour shell(s) too (pair lists = neighbours within the cutoff).
                    for (int sh = 0; sh < config.swap_inject_freeze_shell; ++sh) {
                        std::vector<std::uint8_t> next = inj_free_mask;
                        for (std::size_t k = 0; k < N; ++k) {
                            if (!inj_free_mask[k]) continue;
                            const std::uint32_t cnt = pairs_data[pair_offsets[k]];
                            for (std::uint32_t q = 1; q <= cnt; ++q) next[pairs_data[pair_offsets[k] + q]] = 1;
                        }
                        inj_free_mask.swap(next);
                    }
                    inj_freeze_until = step + static_cast<std::size_t>(config.swap_inject_freeze);
                    // trial state: the pairs of THIS step (for a revert) and the participants' pre-exchange
                    // positions/moments (the injection changed diameters only, so x/m/v here are pre-state);
                    // U was computed by this step's force pass BEFORE the injection -> U0.
                    inj_pairs_saved.clear(); inj_participants.clear();
                    for (std::size_t q = 0; q + 2 < config.swap_inject.size(); q += 3)
                        if (static_cast<std::size_t>(config.swap_inject[q]) == step)
                            inj_pairs_saved.emplace_back(static_cast<std::uint32_t>(config.swap_inject[q + 1]),
                                                         static_cast<std::uint32_t>(config.swap_inject[q + 2]));
                    for (std::size_t k = 0; k < N; ++k) if (inj_free_mask[k]) inj_participants.push_back(static_cast<std::uint32_t>(k));
                    inj_pre_x.clear(); inj_pre_m.clear(); inj_pre_v.clear(); inj_pre_D0.clear();
                    for (std::uint32_t k : inj_participants) inj_pre_D0.push_back(inj_D0_before[k]);
                    inj_pre_sc.clear(); for (std::uint32_t k : inj_participants) inj_pre_sc.push_back(swap_engine.swap_count(k));
                    inj_pre_geo = swap_engine.geo_trades_; inj_pre_cyc = swap_engine.geo_cycle_trades_; inj_pre_tot = swap_total_accepted;
                    inj_free_origin.assign(N, 0);
                    for (std::size_t k = 0; k < N; ++k) if (inj_free_mask[k]) inj_free_origin[particle_origin[k]] = 1;
                    inj_part_origin.clear(); for (std::uint32_t k : inj_participants) inj_part_origin.push_back(static_cast<std::uint32_t>(particle_origin[k]));
                    inj_pairs_origin.clear(); for (const auto& pr : inj_pairs_saved) inj_pairs_origin.emplace_back(static_cast<std::uint32_t>(particle_origin[pr.first]), static_cast<std::uint32_t>(particle_origin[pr.second]));
                    for (std::uint32_t k : inj_participants)
                        for (std::size_t d = 0; d < Ndim; ++d) {
                            inj_pre_x.push_back(x[k * Ndim + d]); inj_pre_m.push_back(m_flat[k * Ndim + d]); inj_pre_v.push_back(v_flat[k * Ndim + d]);
                        }
                    inj_U0 = U; ++inj_trials;
                }
            }
        }
        // Cycle 20: missed-pair sentinel — validates the candidate lists the
        // forces above just consumed. Off unless RCP_SENTINEL is set.
        if (rcp_sentinel_every() > 0 &&
            step % static_cast<std::size_t>(rcp_sentinel_every()) == 0) {
            rcp_sentinel_check(N, Ndim, x.data(), D, box, walls,
                               max_neighbors_per_particle.data(),
                               pair_offsets.data(), pairs_data.data(), step,
                               x_last.data(), D_last.data(), refresh.data(),
                               Dmin * 0.25 * trigger_scale_run);
        }
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
        // Lightweight multi-axis progress log (opt-in via RCP_PROGRESS_STRIDE>0,
        // parsed once). stderr only, no effect on physics/history -> swap=0 stays
        // byte-identical when off. Lets a watcher tail step/phi/overlap/trades to
        // tell legitimate slow convergence from wasted churn [user 2026-09-09].
        static const long g_prog_stride = [](){
            const char* s = std::getenv("RCP_PROGRESS_STRIDE");
            return (s && s[0]) ? std::atol(s) : 0L;
        }();
        if (g_prog_stride > 0 &&
            (step % static_cast<std::size_t>(g_prog_stride) == 0 || step == 1)) {
            std::cerr << "[prog] step=" << step
                      << " phi=" << (phi / phi_modifier)
                      << " maxov=" << max_min_dist
                      << " muflag=" << mu_flag
                      << " geo=" << swap_engine.geo_trades_
                      << " cyc=" << swap_engine.geo_cycle_trades_
                      << " cyc4=" << swap_engine.geo_cycle4_trades_
                      << " capk=" << swap_engine.geo_capk_trades_
                      << " cage=" << swap_engine.geo_cage_trades_
                      << " bneck=" << swap_engine.geo_bottleneck_trades_
                      << "\n";
        }
        if (config.record_visibility && step - 1 < zpos_hist.size()) {
            std::size_t nz = 0;
            for (std::size_t q = 0; q < N; ++q) nz += (z[q] >= 1);
            zpos_hist[step - 1] = static_cast<float>(nz) / static_cast<float>(N);
            mueff_hist[step - 1] = static_cast<float>(mu_eff);
        }
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
        // [207] end of a quench: this step's force pass saw the quenched configuration -> decide.
        if (!inj_free_mask.empty()) inj_remap_mask();                                      // [207] survive reorders
        if (!inj_free_mask.empty() && step == inj_freeze_until + 1 && config.swap_inject_accept) {
            const bool keep = U <= inj_U0 * (1.0 + config.swap_inject_tol);
            std::cerr << "[inj] step " << step << " trial: U0=" << inj_U0 << " U=" << U
                      << (keep ? " ACCEPT" : " REVERT") << "\n";
            if (!keep) {
                ++inj_reverts;
                static thread_local std::vector<std::uint32_t> rv_touched; rv_touched.clear();
                inj_idx_of_origin.assign(N, 0);
                for (std::size_t k = 0; k < N; ++k) inj_idx_of_origin[particle_origin[k]] = static_cast<std::uint32_t>(k);
                for (const auto& pr : inj_pairs_origin)
                    swap_engine.inject(N, D, D0, rv_touched, inj_idx_of_origin[pr.first], inj_idx_of_origin[pr.second]);
                std::size_t q = 0, n_d0_mismatch = 0, n_traded = 0; double rmin = 1e300, rmax = -1e300;
                for (std::size_t pi = 0; pi < inj_part_origin.size(); ++pi) {
                    const std::uint32_t k = inj_idx_of_origin[inj_part_origin[pi]];
                    if (swap_engine.swap_count(k) != inj_pre_sc[pi]) ++n_traded;
                    const double ratio = D0[k] / inj_pre_D0[pi];
                    rmin = std::min(rmin, ratio); rmax = std::max(rmax, ratio);
                    if (std::fabs(ratio - 1.0) > 1e-12) ++n_d0_mismatch;
                    for (std::size_t d = 0; d < Ndim; ++d, ++q) {
                        x[k * Ndim + d] = inj_pre_x[q]; m_flat[k * Ndim + d] = inj_pre_m[q]; v_flat[k * Ndim + d] = inj_pre_v[q];
                    }
                }
                std::cerr << "[inj]   revert check: participants=" << inj_participants.size()
                          << " D0 mismatches after swap-back=" << n_d0_mismatch
                          << " traded-since-arming=" << n_traded
                          << " dgeo=" << (swap_engine.geo_trades_ - inj_pre_geo) << " dcyc=" << (swap_engine.geo_cycle_trades_ - inj_pre_cyc)
                          << " dtotal=" << (swap_total_accepted - inj_pre_tot)
                          << " ratio range [" << rmin << ", " << rmax << "]\n";
                ++g_rcp_perm_epoch;
                swap_need_relayout = true;
            }
            inj_free_mask.clear();
        }
        const bool inj_quench = (step <= inj_freeze_until) && !inj_free_mask.empty();   // [207]
        if (inj_quench) { frz_x = x; frz_m = m_flat; frz_v = v_flat; frz_kappa = kappa; frz_mk = m_kappa; frz_vk = v_kappa; }
        g_t_adam.begin();
        if (rcp_fused_adam() && method != "Verlet") {
            // Cycle 21c: particle ADAM is fused into the position pass
            // below (identical per-element math, m_hat/v_hat carried in
            // registers instead of round-tripping through memory). Only
            // the kappa half of adam_update runs here, because the kappa
            // position update between this site and the fused pass
            // consumes m_hat_kappa / v_hat_kappa. Lines are verbatim from
            // adam_update.
            m_kappa = beta1 * m_kappa - (1 - beta1) * dkappa;
            v_kappa = beta2 * v_kappa + (1 - beta2) * (dkappa * dkappa);
            v_update_kappa = v_kappa;
            m_hat_kappa = m_kappa / (1.0 - std::pow(beta1, double(t)));
            v_hat_kappa = v_update_kappa / (1.0 - std::pow(beta2, double(t)));
        } else {
            adam_update(method, N, Ndim,
                        F, dkappa, beta1, beta2, t, dt, verlet_drag,
                        m_flat.data(), v_flat.data(), v_max,
                        v_update_flat.data(),
                        m_hat_flat.data(), v_hat_flat.data(),
                        a, v_verlet, a_old,
                        m_kappa, v_kappa, v_update_kappa,
                        m_hat_kappa, v_hat_kappa);
        }
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

        if (!runtime_fix_diameter && !quench_on) {   // [221] kappa frozen during a quench
            // [219] while a morph is advancing, damp the global-scale response so positions absorb the transient
            // imbalance of the reference-size shift (morph_kappa_damp = 1 -> original step, bit for bit)
            const double kd = (morph_on && morph_armed && morph_lambda < 1.0) ? config.morph_kappa_damp : 1.0;
            kappa = kappa - kd * alpha * m_hat_kappa / (std::sqrt(v_hat_kappa) + epsilon);   // kd == 1.0 is exact
        }
        g_t_bookkeep.end();
        g_t_position.begin();
        if (quench_on && ell_on) {
            // ---- ellipsoid quench (additive). The same FIRE minimisation as the sphere branch below, on positions
            // AND orientations. A particle's generalised force is (F_i, T_i / a_i) and its generalised velocity
            // (v_i, a_i w_i), with a_i the bounding radius: every component is then a length (or a force), and the
            // power F.v + T.w is unchanged. Residual of a load-bearing particle: |(F_i, T_i/a_i)| / S_i, S_i the sum
            // of its contact force magnitudes, so it measures force AND torque balance. Not implemented here: the
            // two-level, activity-damping, warm-up and per-particle-reset variants (refused at input).
            const std::size_t nrot = ell_nrot(Ndim), nn = Ndim * Ndim, nd = N * Ndim;
            if (q_iter == 0) {
                std::fill(q_held.begin(), q_held.end(), 0);
                if (!config.quench_held.empty()) {
                    std::vector<char> want(N, 0);
                    for (double hi : config.quench_held) { const long k = static_cast<long>(hi); if (k >= 0 && k < static_cast<long>(N)) want[k] = 1; }
                    for (std::size_t k = 0; k < N; ++k) q_held[k] = want[particle_origin[k]];
                }
                q_any_held = std::any_of(q_held.begin(), q_held.end(), [](char h) { return h != 0; });
                if (config.quench_kappa_mult != 1.0 && config.quench_kappa_mult > 0.0) kappa *= config.quench_kappa_mult;
            }
            // flat-wall contacts: half-width h = sqrt(A_dd) in the role of the sphere's radius, force 2 (h - gap)
            const bool q_walls = std::any_of(walls.begin(), walls.end(), [](std::int8_t w) { return w != 0; });
            if (q_walls) {
                q_wall_s.assign(N, 0.0); q_wall_n.assign(N, 0);
                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t d = 0; d < Ndim; ++d) {
                        if (walls[d] != 1) continue;
                        const double h = std::sqrt(ell_amat[i * nn + d * Ndim + d]);
                        for (int side = 0; side < 2; ++side) {
                            const double gap = std::abs(box[d] * side - x[i * Ndim + d]);
                            if (gap < h) { q_wall_s[i] += 2.0 * (h - gap); ++q_wall_n[i]; }
                        }
                    }
                    if (walls[0] < 0) {   // curved container: force 2 (h - delta) at the farthest point, one contact
                        CurvedWallSum cw;
                        const double Rc = 0.5 * box[0];
                        if (shape_curved_wall(ell_shape, Ndim, curved_axes(walls[0], Ndim), Rc, 1.0, x.data() + i * Ndim, D[i],
                                              ell_amat.data() + i * nn, ell_rot.data() + i * nn, ell_q.data() + i * Ndim,
                                              nullptr, nullptr, cw)) { q_wall_s[i] += cw.f; q_wall_n[i] += cw.count; }
                    }
                }
            }
            auto qS = [&](std::size_t i) { return q_walls ? swap_sink.gsum[i] + q_wall_s[i] : swap_sink.gsum[i]; };
            auto qZ = [&](std::size_t i) { return q_walls ? static_cast<long>(z[i]) + q_wall_n[i] : static_cast<long>(z[i]); };
            auto qA = [&](std::size_t i) { return 0.5 * D[i]; };                          // bounding radius
            const long q_zmin = static_cast<long>(Ndim) + 1;
            auto gen_force2 = [&](std::size_t i) {
                double g2 = 0.0;
                for (std::size_t d = 0; d < Ndim; ++d) g2 += F[i * Ndim + d] * F[i * Ndim + d];
                const double ia = 1.0 / qA(i);
                for (std::size_t c = 0; c < nrot; ++c) { const double g = ell_torque[i * nrot + c] * ia; g2 += g * g; }
                return g2;
            };
            double f2 = 0.0, ssum = 0.0; long zsum = 0;
            for (std::size_t d = 0; d < Ndim; ++d) q_held_force[d] = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) { for (std::size_t d = 0; d < Ndim; ++d) q_held_force[d] += F[i * Ndim + d]; continue; }
                const double S = qS(i);
                if (S > 0.0) { ssum += S; zsum += qZ(i); }
            }
            const double fmean = zsum > 0 ? ssum / static_cast<double>(zsum) : 1.0;
            const double s_floor = 1e-2 * fmean * 4.0;
            double rmax = 0.0; std::size_t nback = 0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) continue;
                const double S = qS(i);
                if (qZ(i) >= q_zmin && S >= s_floor) {
                    const double g2 = gen_force2(i);
                    ++nback; f2 += g2; const double r = std::sqrt(g2) / S; if (r > rmax) rmax = r;
                }
            }
            q_res = rmax; q_gres = (fmean > 0.0 && nback > 0) ? std::sqrt(f2 / static_cast<double>(std::max<std::size_t>(nback, 1))) / fmean : 0.0;
            if (q_iter % std::max<long>(1, config.quench_report_every) == 0) {
                q_log.push_back(static_cast<double>(q_iter)); q_log.push_back(q_res); q_log.push_back(q_gres);
                q_log.push_back(U); q_log.push_back(static_cast<double>(nback)); q_log.push_back(q_dt);
            }
            const bool q_armed = (q_iter >= 1) || (config.quench_kappa_mult == 1.0);
            if (q_armed && nback > 0 && q_res <= config.quench_tol) q_status = 1;
            else if (!std::isfinite(q_res) || !std::isfinite(U) || !std::isfinite(q_gres)) q_status = 6;
            else if (q_armed && nback == 0 && (U == 0.0 || ++q_unjam_streak >= 50)) q_status = 4;
            else if (q_iter >= config.quench_max_iter) q_status = 2;
            if (nback > 0) q_unjam_streak = 0;
            if (q_status == 0 && q_armed && nback > 0) {
                // progress = the residual halved (as for spheres) OR the energy fell by more than 0.1 % since the
                // reference point: ellipsoid packings rearrange slowly, lowering U for tens of thousands of
                // iterations while the largest per-particle residual wanders
                if (q_stall_ref < 0.0 || q_res <= config.quench_stall_factor * q_stall_ref || U < 0.999 * ell_q_stall_u) { q_stall_ref = q_res; q_stall_ref_iter = q_iter; ell_q_stall_u = U; }
                else if (config.quench_stall_window > 0 && q_iter - q_stall_ref_iter >= config.quench_stall_window) q_status = 3;
                if (q_status == 0 && config.quench_dt_min > 0.0 && q_dt < config.quench_dt_min) q_status = 5;
            }
            if (q_status != 0) {
                q_converged = (q_status == 1); steps = step; g_t_position.end(); break;
            }
            if (config.quench_mass_scaled && (q_iter % 200 == 0 || q_minv.size() != N)) {
                q_minv.assign(N, 1.0); std::vector<double> zz; zz.reserve(N);
                for (std::size_t i = 0; i < N; ++i) if (qZ(i) > 0) zz.push_back(static_cast<double>(qZ(i)));
                double med = 1.0; if (!zz.empty()) { std::nth_element(zz.begin(), zz.begin() + zz.size() / 2, zz.end()); med = zz[zz.size() / 2]; }
                for (std::size_t i = 0; i < N; ++i) q_minv[i] = med / std::max(1.0, static_cast<double>(qZ(i)));
            }
            auto q_free = [&](std::size_t i) { return !(q_held[i] || qS(i) <= 0.0); };
            // FIRE: P = G.V over translational and rotational components
            double P = 0.0, fnorm2 = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                if (!q_free(i)) continue;
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0, ia = 1.0 / qA(i);
                for (std::size_t d = 0; d < Ndim; ++d) { const double a = F[i * Ndim + d] * mi; P += a * q_v[i * Ndim + d]; fnorm2 += a * a; }
                for (std::size_t c = 0; c < nrot; ++c) { const double a = ell_torque[i * nrot + c] * ia * mi; P += a * ell_qw[i * nrot + c]; fnorm2 += a * a; }
            }
            if (P > 0.0) { if (++q_npos > 5) { q_dt = std::min(q_dt * 1.1, config.quench_dt_max); q_alpha *= 0.99; } }
            else {
                q_npos = 0; q_dt *= 0.5; q_alpha = config.quench_alpha0;
                // half a step back, in position and in orientation, then stop
                for (std::size_t k = 0; k < nd; ++k) { x[k] -= 0.5 * q_dt * q_v[k]; q_v[k] = 0.0; }
                double wback[8 * 7 / 2];
                for (std::size_t i = 0; i < N; ++i) {
                    const double ia = 1.0 / qA(i); bool any = false;
                    for (std::size_t c = 0; c < nrot; ++c) { wback[c] = -0.5 * q_dt * ell_qw[i * nrot + c] * ia; if (wback[c] != 0.0) any = true; ell_qw[i * nrot + c] = 0.0; }
                    if (any) ell_rotate_by(Ndim, ell_rot.data() + i * nn, wback);
                }
            }
            const double fn = std::sqrt(fnorm2);
            // pass 1: V <- V + dt G/m
            double vnew2 = 0.0; std::size_t nact = 0;
            for (std::size_t i = 0; i < N; ++i) {
                if (!q_free(i)) {
                    for (std::size_t d = 0; d < Ndim; ++d) q_v[i * Ndim + d] = 0.0;
                    for (std::size_t c = 0; c < nrot; ++c) ell_qw[i * nrot + c] = 0.0;
                    continue;
                }
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0, ia = 1.0 / qA(i); ++nact;
                for (std::size_t d = 0; d < Ndim; ++d) { const double vv = q_v[i * Ndim + d] + q_dt * F[i * Ndim + d] * mi; q_v[i * Ndim + d] = vv; vnew2 += vv * vv; }
                for (std::size_t c = 0; c < nrot; ++c) { const double vv = ell_qw[i * nrot + c] + q_dt * ell_torque[i * nrot + c] * ia * mi; ell_qw[i * nrot + c] = vv; vnew2 += vv * vv; }
            }
            const double vn = std::sqrt(vnew2);
            // pass 2: mixing V <- (1 - alpha) V + alpha |V| G/|G| (global norms, or per particle for mix_mode 1)
            std::vector<double> vmean(Ndim, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                if (!q_free(i)) continue;
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0, ia = 1.0 / qA(i);
                double vn_i = vn, fn_i = fn;
                if (config.quench_mix_mode >= 1) {
                    double a2 = 0.0, b2 = 0.0;
                    for (std::size_t d = 0; d < Ndim; ++d) { const double a = F[i * Ndim + d] * mi; a2 += a * a; b2 += q_v[i * Ndim + d] * q_v[i * Ndim + d]; }
                    for (std::size_t c = 0; c < nrot; ++c) { const double a = ell_torque[i * nrot + c] * ia * mi; a2 += a * a; b2 += ell_qw[i * nrot + c] * ell_qw[i * nrot + c]; }
                    fn_i = std::sqrt(a2); vn_i = std::sqrt(b2);
                }
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double vv = q_v[i * Ndim + d];
                    if (fn_i > 0.0) vv = (1.0 - q_alpha) * vv + q_alpha * vn_i * (F[i * Ndim + d] * mi) / fn_i;
                    q_v[i * Ndim + d] = vv; vmean[d] += vv;
                }
                for (std::size_t c = 0; c < nrot; ++c) {
                    double vv = ell_qw[i * nrot + c];
                    if (fn_i > 0.0) vv = (1.0 - q_alpha) * vv + q_alpha * vn_i * (ell_torque[i * nrot + c] * ia * mi) / fn_i;
                    ell_qw[i * nrot + c] = vv;
                }
            }
            for (std::size_t d = 0; d < Ndim; ++d) vmean[d] = nact ? vmean[d] / static_cast<double>(nact) : 0.0;
            // move: positions (uniform translation removed along periodic directions when nothing is held), then orientations
            double wstep[8 * 7 / 2];
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) continue;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    if (qS(i) > 0.0 && walls[d] == 0 && !q_any_held) q_v[i * Ndim + d] -= vmean[d];
                    double xv = x[i * Ndim + d] + q_dt * q_v[i * Ndim + d]; const double L = box[d];
                    if (xv >= L) xv -= L; else if (xv < 0.0) xv += L; x[i * Ndim + d] = xv;
                }
                const double ia = 1.0 / qA(i); bool any = false;
                for (std::size_t c = 0; c < nrot; ++c) { wstep[c] = q_dt * ell_qw[i * nrot + c] * ia; if (wstep[c] != 0.0) any = true; }
                if (any) ell_rotate_by(Ndim, ell_rot.data() + i * nn, wstep);
            }
            ++q_iter;
        } else
        if (quench_on) {   // [221] FIRE quench: adaptive inertial descent on the true energy, unit or z-scaled masses
            const std::size_t nd = N * Ndim;
            if (q_iter == 0) {
                // arm: the held mask refers to the slot order the caller sees (the resumed blob). It must be built HERE,
                // after the resume, because the initialiser spatially reorders the fresh arrays before the checkpoint
                // state overwrites them ([222]: a mask built at declaration ended up permuted against the particles).
                // quench_held holds ORIGIN ids (indices into the arrays the caller handed in; for a checkpoint chain,
                // the first run's order): particle_origin[k] is the stable identity that survives every spatial
                // reorder and every resume, so the mask is built through it and then permuted with the slots.
                std::fill(q_held.begin(), q_held.end(), 0);
                if (!config.quench_held.empty()) {
                    std::vector<char> want(N, 0);
                    for (double hi : config.quench_held) { const long k = static_cast<long>(hi); if (k >= 0 && k < static_cast<long>(N)) want[k] = 1; }
                    for (std::size_t k = 0; k < N; ++k) q_held[k] = want[particle_origin[k]];
                }
                q_any_held = std::any_of(q_held.begin(), q_held.end(), [](char h) { return h != 0; });
                if (config.quench_kappa_mult != 1.0 && config.quench_kappa_mult > 0.0) kappa *= config.quench_kappa_mult;   // takes effect at the next D = D0*kappa
            }
            // [225] wall contacts (planar walls at 0 and L in each dimension with walls[d] == 1; a curved container on
            // the first -walls[0] dimensions): the force kernel's mirror-image law gives |f| = 2 (r - gap distance).
            // They count toward the contact number and the contact-force sum exactly as a pair contact does.
            const bool q_walls = std::any_of(walls.begin(), walls.end(), [](std::int8_t w) { return w != 0; });
            if (q_walls) {
                q_wall_s.assign(N, 0.0); q_wall_n.assign(N, 0);
                const std::size_t t_curved = walls[0] < 0 ? static_cast<std::size_t>(-walls[0]) : 0;
                for (std::size_t i = 0; i < N; ++i) {
                    const double r = D[i] / 2.0;
                    for (std::size_t d = 0; d < Ndim; ++d) {
                        if (walls[d] != 1) continue;
                        for (int side = 0; side < 2; ++side) {
                            const double gap = std::abs(box[d] * side - x[i * Ndim + d]);
                            if (gap < r) { q_wall_s[i] += 2.0 * (r - gap); ++q_wall_n[i]; }
                        }
                    }
                    if (t_curved > 0) {
                        const double R = box[0] / 2.0; double dc2 = 0.0;
                        for (std::size_t m = 0; m < t_curved && m < Ndim; ++m) { const double u = x[i * Ndim + m] - R; dc2 += u * u; }
                        const double gap = R - std::sqrt(dc2);   // the kernel's test; its force uses |gap|
                        if (gap < r) { q_wall_s[i] += 2.0 * std::abs(r - std::abs(gap)); ++q_wall_n[i]; }
                    }
                }
            }
            // contact-force sum S_i (K = 1: overlap sum = force sum) and contact number, walls included
            auto qS = [&](std::size_t i) { return q_walls ? swap_sink.gsum[i] + q_wall_s[i] : swap_sink.gsum[i]; };
            auto qZ = [&](std::size_t i) { return q_walls ? static_cast<long>(z[i]) + q_wall_n[i] : static_cast<long>(z[i]); };
            const long q_zmin = static_cast<long>(Ndim) + 1;   // a particle needs Ndim + 1 contacts to be caged (4 in 3D)
            // residuals on the backbone from F and the per-particle contact-force sums
            // load-bearing backbone: z >= Ndim + 1 AND contact-force sum >= 1e-2 x (mean contact force x 4): below that
            // the ratio |F_i|/S_i is limited by position rounding (K*1e-16 / S_i), not by balance ([222])
            double f2 = 0.0, ssum = 0.0; long zsum = 0;
            for (std::size_t d = 0; d < Ndim; ++d) q_held_force[d] = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) { for (std::size_t d = 0; d < Ndim; ++d) q_held_force[d] += F[i * Ndim + d]; continue; }
                const double S = qS(i);
                if (S > 0.0) { ssum += S; zsum += qZ(i); }
            }
            const double fmean = zsum > 0 ? ssum / static_cast<double>(zsum) : 1.0;
            const double s_floor = 1e-2 * fmean * 4.0;
            double rmax = 0.0; std::size_t nback = 0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) continue;
                const double S = qS(i);
                if (qZ(i) >= q_zmin && S >= s_floor) {
                    double fn2 = 0.0; for (std::size_t d = 0; d < Ndim; ++d) fn2 += F[i * Ndim + d] * F[i * Ndim + d];
                    ++nback; f2 += fn2; const double r = std::sqrt(fn2) / S; if (r > rmax) rmax = r;
                }
            }
            q_res = rmax; q_gres = (fmean > 0.0 && nback > 0) ? std::sqrt(f2 / static_cast<double>(std::max<std::size_t>(nback, 1))) / fmean : 0.0;
            if (q_iter % std::max<long>(1, config.quench_report_every) == 0) {
                q_log.push_back(static_cast<double>(q_iter)); q_log.push_back(q_res); q_log.push_back(q_gres);
                q_log.push_back(U); q_log.push_back(static_cast<double>(nback)); q_log.push_back(q_dt);
            }
            // a kappa multiplier takes effect at the NEXT D = D0*kappa: never declare convergence before one iteration
            const bool q_armed = (q_iter >= 1) || (config.quench_kappa_mult == 1.0);
            // [225] stopping rules, first match wins: converged (needs a non-empty backbone - an empty one is vacuous),
            // non-finite, unjammed (no overlap at all, or an empty backbone for 50 consecutive iterations), iteration
            // cap, stall (best residual not reduced by stall_factor within stall_window iterations), step collapse
            // (FIRE dt below dt_min while a backbone exists: FIRE finds no step it can follow downhill)
            if (q_armed && nback > 0 && q_res <= config.quench_tol) q_status = 1;
            else if (!std::isfinite(q_res) || !std::isfinite(U) || !std::isfinite(q_gres)) q_status = 6;
            else if (q_armed && nback == 0 && (U == 0.0 || ++q_unjam_streak >= 50)) q_status = 4;
            else if (q_iter >= config.quench_max_iter) q_status = 2;
            if (nback > 0) q_unjam_streak = 0;
            if (q_status == 0 && q_armed && nback > 0) {
                if (q_stall_ref < 0.0 || q_res <= config.quench_stall_factor * q_stall_ref) { q_stall_ref = q_res; q_stall_ref_iter = q_iter; }
                else if (config.quench_stall_window > 0 && q_iter - q_stall_ref_iter >= config.quench_stall_window) q_status = 3;
                if (q_status == 0 && config.quench_dt_min > 0.0 && q_dt < config.quench_dt_min) q_status = 5;
            }
            if (q_status != 0) {
                q_converged = (q_status == 1); steps = step; g_t_position.end(); break;
            }
            // [223] two-level classification by load: heavy = caged (z >= Ndim + 1) AND load >= light_frac x typical; refreshed
            // every class_every iterations (not every step: a flickering set would perturb FIRE's power sum)
            if (config.quench_two_level && (q_iter % std::max<int>(1, config.quench_class_every) == 0)) {
                const double s_light = config.quench_light_frac * fmean * 4.0;
                for (std::size_t i = 0; i < N; ++i) q_light[i] = (q_held[i] || qZ(i) < q_zmin || qS(i) < s_light) ? 1 : 0;
            } else if (!config.quench_two_level) { std::fill(q_light.begin(), q_light.end(), 0); }
            // masses
            if (config.quench_mass_scaled && (q_iter % 200 == 0 || q_minv.size() != N)) {
                q_minv.assign(N, 1.0); std::vector<double> zz; zz.reserve(N);
                for (std::size_t i = 0; i < N; ++i) if (qZ(i) > 0) zz.push_back(static_cast<double>(qZ(i)));
                double med = 1.0; if (!zz.empty()) { std::nth_element(zz.begin(), zz.begin() + zz.size() / 2, zz.end()); med = zz[zz.size() / 2]; }
                for (std::size_t i = 0; i < N; ++i) q_minv[i] = med / std::max(1.0, static_cast<double>(qZ(i)));
            }
            // FIRE: P = F.v ; mixing ; adaptive dt / alpha ; semi-implicit Euler ; held and translation handling
            double P = 0.0, vnorm2 = 0.0, fnorm2 = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i] || q_light[i] || qS(i) <= 0.0) continue;
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0;
                for (std::size_t d = 0; d < Ndim; ++d) { const double a = F[i * Ndim + d] * mi; P += a * q_v[i * Ndim + d]; fnorm2 += a * a; vnorm2 += q_v[i * Ndim + d] * q_v[i * Ndim + d]; }
            }
            const bool q_warm = q_iter < config.quench_warm_steps;   // [224] warm phase: no resets, isokinetic rescale below
            if (P > 0.0 || q_warm) { if (++q_npos > 5) { q_dt = std::min(q_dt * 1.1, config.quench_dt_max); q_alpha *= 0.99; } }
            else { q_npos = 0; q_dt *= 0.5; q_alpha = config.quench_alpha0;
                   for (std::size_t k = 0; k < nd; ++k) { x[k] -= 0.5 * q_dt * q_v[k]; q_v[k] = 0.0; } vnorm2 = 0.0; }
            if (config.quench_mix_mode >= 2) {
                // per-particle reset [222]: a particle moving against its own force stops (half-step back), independent of the rest
                for (std::size_t i = 0; i < N; ++i) {
                    if (q_held[i] || q_light[i] || qS(i) <= 0.0) continue;
                    double Pi = 0.0; for (std::size_t d = 0; d < Ndim; ++d) Pi += F[i * Ndim + d] * q_v[i * Ndim + d];
                    if (Pi < 0.0) {
                        if (config.quench_mix_mode == 3) { for (std::size_t d = 0; d < Ndim; ++d) q_v[i * Ndim + d] *= 0.5; }   // [224] soft per-particle damping: local dissipation without the discontinuity of a reset
                        else for (std::size_t d = 0; d < Ndim; ++d) { x[i * Ndim + d] -= 0.5 * q_dt * q_v[i * Ndim + d]; q_v[i * Ndim + d] = 0.0; }
                    }
                }
            }
            const double fn = std::sqrt(fnorm2);
            // pass 1: v <- v + dt a (held / force-free particles carry no velocity); accumulate |v| after the update
            double vnew2 = 0.0; std::size_t nact = 0;
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i] || q_light[i] || qS(i) <= 0.0) { for (std::size_t d = 0; d < Ndim; ++d) q_v[i * Ndim + d] = 0.0; continue; }
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0; ++nact;
                for (std::size_t d = 0; d < Ndim; ++d) { const double vv = q_v[i * Ndim + d] + q_dt * F[i * Ndim + d] * mi; q_v[i * Ndim + d] = vv; vnew2 += vv * vv; }
            }
            const double vn = std::sqrt(vnew2);
            // pass 2: FIRE mixing v <- (1-alpha) v + alpha |v| F/|F|, then the translation mean
            // [228] one entry per dimension (was double[3]: out of bounds for Ndim >= 4 -> collective drift)
            std::vector<double> vmean(Ndim, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i] || q_light[i] || qS(i) <= 0.0) continue;
                const double mi = config.quench_mass_scaled ? q_minv[i] : 1.0;
                double vn_i = vn, fn_i = fn;
                if (config.quench_mix_mode >= 1) {   // per-particle norms [222]: no global velocity scale is imposed on a lightly loaded particle
                    double a2 = 0.0, b2 = 0.0;
                    for (std::size_t d = 0; d < Ndim; ++d) { const double a = F[i * Ndim + d] * mi; a2 += a * a; b2 += q_v[i * Ndim + d] * q_v[i * Ndim + d]; }
                    fn_i = std::sqrt(a2); vn_i = std::sqrt(b2);
                }
                for (std::size_t d = 0; d < Ndim; ++d) {
                    double vv = q_v[i * Ndim + d];
                    if (fn_i > 0.0) vv = (1.0 - q_alpha) * vv + q_alpha * vn_i * (F[i * Ndim + d] * mi) / fn_i;
                    q_v[i * Ndim + d] = vv; vmean[d] += vv;
                }
            }
            for (std::size_t d = 0; d < Ndim; ++d) vmean[d] = nact ? vmean[d] / static_cast<double>(nact) : 0.0;
            if (config.quench_freeze_damp > 0.0) {
                // [224] activity-gated dissipation: EMA of |v| per particle; particles far below the median activity are
                // damped each step, active clusters untouched
                if (q_act.size() != N) q_act.assign(N, 0.0);
                std::vector<double> tmp; tmp.reserve(N);
                for (std::size_t i = 0; i < N; ++i) {
                    double v2 = 0.0; for (std::size_t d = 0; d < Ndim; ++d) v2 += q_v[i * Ndim + d] * q_v[i * Ndim + d];
                    q_act[i] = 0.99 * q_act[i] + 0.01 * std::sqrt(v2);
                    if (!(q_held[i] || q_light[i] || qS(i) <= 0.0)) tmp.push_back(q_act[i]);
                }
                if (q_iter % 200 == 0 && !tmp.empty()) { std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end()); q_act_med = tmp[tmp.size() / 2]; }
                if (q_act_med > 0.0 && q_iter > 400) {
                    const double thr = config.quench_freeze_frac * q_act_med, keep = 1.0 - config.quench_freeze_damp;
                    for (std::size_t i = 0; i < N; ++i) if (q_act[i] < thr) for (std::size_t d = 0; d < Ndim; ++d) q_v[i * Ndim + d] *= keep;
                }
            }
            if (q_warm && nact > 0 && fmean > 0.0) {
                // [224] isokinetic thermostat: mean speed such that a particle moves quench_warm_amp mean overlaps per step
                double v2 = 0.0;
                for (std::size_t i = 0; i < N; ++i) { if (q_held[i] || q_light[i] || qS(i) <= 0.0) continue; for (std::size_t d = 0; d < Ndim; ++d) v2 += q_v[i * Ndim + d] * q_v[i * Ndim + d]; }
                const double vrms = std::sqrt(v2 / static_cast<double>(nact)); const double vt = config.quench_warm_amp * fmean / std::max(q_dt, 1e-300);
                if (vrms > 0.0) { const double sc = vt / vrms; for (std::size_t k = 0; k < nd; ++k) q_v[k] *= sc; for (std::size_t d = 0; d < Ndim; ++d) vmean[d] *= sc; }
            }
            for (std::size_t i = 0; i < N; ++i) {
                if (q_held[i]) continue;
                if (q_light[i]) {
                    // [223] light particle: overdamped local step against the (frozen this step) heavy set; F/(K z) is the exact
                    // opening displacement for one contact and an under-relaxed Jacobi step for a cage; no velocity, no FIRE state
                    const double zi = static_cast<double>(qZ(i));
                    if (zi <= 0.0 || qS(i) <= 0.0) continue;
                    const double fac = config.quench_light_omega / zi;   // K = 1 in this kernel (overlap sum = contact-force sum)
                    for (std::size_t d = 0; d < Ndim; ++d) {
                        double xv = x[i * Ndim + d] + fac * F[i * Ndim + d]; const double L = box[d];
                        if (xv >= L) xv -= L; else if (xv < 0.0) xv += L; x[i * Ndim + d] = xv;
                    }
                    continue;
                }
                for (std::size_t d = 0; d < Ndim; ++d) {
                    // [225] the uniform translation is a zero mode only along periodic directions and only when nothing
                    // is held: walls and held particles exert a physical net force, and removing the response to it
                    // leaves every free particle a common residual force (FIRE then stalls with its step collapsed)
                    if (qS(i) > 0.0 && walls[d] == 0 && !q_any_held) q_v[i * Ndim + d] -= vmean[d];
                    double xv = x[i * Ndim + d] + q_dt * q_v[i * Ndim + d]; const double L = box[d];
                    if (xv >= L) xv -= L; else if (xv < 0.0) xv += L; x[i * Ndim + d] = xv;
                }
            }
            ++q_iter;
        } else {
        if (rcp_fused_adam() && method != "Verlet") {
            // Cycle 21c: fused particle-ADAM + position update + PBC wrap.
            // Per-element math is verbatim adam_update (m/v EMA + bias
            // correction) followed by the position step and the Cycle-15k
            // branch wrap — one pass over F/m/v/x instead of three passes
            // with m_hat/v_hat round-tripped through memory. The sign-flip
            // pre-pass (when enabled) already ran before this point, and
            // nothing between the original adam site and here touches
            // m_flat/v_flat, so ordering semantics are unchanged.
            const double inv_b1_corr = 1.0 / (1.0 - std::pow(beta1, double(t)));
            const double inv_b2_corr = 1.0 / (1.0 - std::pow(beta2, double(t)));
            const double* Fd = F;
            double* md = m_flat.data();
            double* vd = v_flat.data();
            double* xd = x.data();
            const double* box_data = box.data();
            const bool use_pp_mult = pp_alpha_mult_loaded && (mu_flag == -1);
            const double* ppm = use_pp_mult ? pp_alpha_mult.data() : nullptr;
#if defined(__AVX512F__)
            if (!use_pp_mult && (Ndim == 2 || Ndim == 3)) {
                // Cycle 21g: 8-wide SIMD over the flat array. Lane math
                // mirrors the scalar fused loop op-for-op (same contraction
                // shapes; vector sqrt/div are correctly-rounded IEEE, so
                // per-lane results are bit-identical), and the branch wrap
                // becomes two masked add/subs. The per-element box length
                // cycles with period Ndim; lanes see rotated L vectors
                // selected by (block offset % Ndim). Parallelized over
                // 8-element blocks, schedule(static) — deterministic
                // assignment, disjoint writes.
                const std::size_t total = N * Ndim;
                const __m512d b1v   = _mm512_set1_pd(beta1);
                const __m512d c1v   = _mm512_set1_pd(1.0 - beta1);
                const __m512d b2v   = _mm512_set1_pd(beta2);
                const __m512d c2v   = _mm512_set1_pd(1.0 - beta2);
                const __m512d mb1c  = _mm512_set1_pd(inv_b1_corr);
                const __m512d mb2c  = _mm512_set1_pd(inv_b2_corr);
                const __m512d alpv  = _mm512_set1_pd(alpha);
                const __m512d epsv  = _mm512_set1_pd(epsilon);
                const __m512d zerov = _mm512_setzero_pd();
                // Rotated box-length vectors: rotL[r][lane] =
                // box[(r + lane) % Ndim].
                alignas(64) double rotL[3][8];
                for (std::size_t r = 0; r < Ndim; ++r)
                    for (int lane = 0; lane < 8; ++lane)
                        rotL[r][lane] =
                            box_data[(r + static_cast<std::size_t>(lane)) % Ndim];
                const std::size_t nblocks = (total + 7) / 8;
                #pragma omp parallel for schedule(static)
                for (std::size_t b = 0; b < nblocks; ++b) {
                    const std::size_t idx = b * 8;
                    const std::size_t remaining = total - idx;
                    const __mmask8 lm =
                        remaining >= 8
                            ? static_cast<__mmask8>(0xFF)
                            : static_cast<__mmask8>((1u << remaining) - 1u);
                    const __m512d Lv =
                        _mm512_load_pd(rotL[idx % Ndim]);
                    const __m512d Fv = _mm512_maskz_loadu_pd(lm, Fd + idx);
                    const __m512d mv = _mm512_maskz_loadu_pd(lm, md + idx);
                    const __m512d vv = _mm512_maskz_loadu_pd(lm, vd + idx);
                    const __m512d xv0 = _mm512_maskz_loadu_pd(lm, xd + idx);
                    // mk = beta1*m - (1-beta1)*F — gcc folds the subtract
                    // into the LEFT product: fmsub(b1, m, c1*F).
                    __m512d mk =
                        _mm512_fmsub_pd(b1v, mv, _mm512_mul_pd(c1v, Fv));
                    // vk = beta2*v + (1-beta2)*F*F — gcc folds the add into
                    // the LEFT product: fmadd(b2, v, (c2*F)*F).
                    __m512d vk = _mm512_fmadd_pd(
                        b2v, vv,
                        _mm512_mul_pd(_mm512_mul_pd(c2v, Fv), Fv));
                    _mm512_mask_storeu_pd(md + idx, lm, mk);
                    _mm512_mask_storeu_pd(vd + idx, lm, vk);
                    const __m512d mh = _mm512_mul_pd(mk, mb1c);
                    const __m512d vh = _mm512_mul_pd(vk, mb2c);
                    const __m512d denom =
                        _mm512_add_pd(_mm512_sqrt_pd(vh), epsv);
                    const __m512d step =
                        _mm512_div_pd(_mm512_mul_pd(alpv, mh), denom);
                    __m512d xv = _mm512_sub_pd(xv0, step);
                    const __mmask8 ge =
                        _mm512_cmp_pd_mask(xv, Lv, _CMP_GE_OQ);
                    xv = _mm512_mask_sub_pd(xv, ge, xv, Lv);
                    const __mmask8 lt0 =
                        _mm512_cmp_pd_mask(xv, zerov, _CMP_LT_OQ);
                    xv = _mm512_mask_add_pd(xv, lt0, xv, Lv);
                    _mm512_mask_storeu_pd(xd + idx, lm, xv);
                }
            } else
#endif
            {
            #pragma omp parallel for schedule(static)
            for (std::size_t k = 0; k < N; ++k) {
                const std::size_t base = k * Ndim;
                const double alpha_k =
                    use_pp_mult ? alpha * ppm[k] : alpha;
                for (std::size_t kk = 0; kk < Ndim; ++kk) {
                    const std::size_t idx = base + kk;
                    const double Fkd = Fd[idx];
                    const double mk = beta1 * md[idx] - (1.0 - beta1) * Fkd;
                    const double vk = beta2 * vd[idx] + (1.0 - beta2) * Fkd * Fkd;
                    md[idx] = mk;
                    vd[idx] = vk;
                    const double mhk = mk * inv_b1_corr;
                    const double vhk = vk * inv_b2_corr;
                    double xv = xd[idx] - alpha_k * mhk / (std::sqrt(vhk) + epsilon);
                    const double L = box_data[kk];
                    if (xv >= L)     xv -= L;
                    else if (xv < 0) xv += L;
                    xd[idx] = xv;
                }
            }
            }
        } else {
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
        }  // end non-fused path
        }  // end quench/adam branch [221]
        if (ell_on && !ell_shape.all_one && !quench_on) {   // ellipsoids: Adam step on the orientations, same alpha, betas and step counter (spheres do not turn)
            ell_rotation_step(N, Ndim, ell_torque.data(), ell_m.data(), ell_v.data(), ell_rot.data(), D.data(),
                              alpha * config.ell_rot_scale, beta1, beta2, double(t), epsilon);
        }
        g_t_position.end();
        if (inj_quench) {                                                                // [207] hold everyone else
            for (std::size_t k = 0; k < N; ++k) {
                if (inj_free_mask[k]) continue;
                for (std::size_t d = 0; d < Ndim; ++d) {
                    x[k * Ndim + d] = frz_x[k * Ndim + d];
                    m_flat[k * Ndim + d] = frz_m[k * Ndim + d];
                    v_flat[k * Ndim + d] = frz_v[k * Ndim + d];
                }
            }
            kappa = frz_kappa; m_kappa = frz_mk; v_kappa = frz_vk;                      // kappa held (size ratios and the energy comparison intact)
        }
        g_t_bookkeep.begin();

        U_history[step - 1] = U;
        hist_filled = step;   // all per-step history series now filled through [0, step)

        const bool should_report =
            observer && progress_interval > 0 && (step % progress_interval == 0);
        // samples every trajectory_interval steps, or (trajectory_per_decade > 0) evenly in log(step)
        bool should_capture =
            capture_positions && trajectory_interval > 0 && (step % trajectory_interval == 0);
        if (capture_positions && config.trajectory_per_decade > 0) {
            should_capture = static_cast<std::size_t>(step) >= traj_log_next;
            if (should_capture) {
                const double grown = static_cast<double>(step) * std::pow(10.0, 1.0 / config.trajectory_per_decade);
                traj_log_next = std::max(static_cast<std::size_t>(step) + 1, static_cast<std::size_t>(std::llround(grown)));
            }
        }
        if (should_capture) traj_last_captured = static_cast<std::size_t>(step);
        if (should_report || should_capture) {
            record_observation(step, should_capture);
        }

        // checkpoint v3: wall-timer write at end-of-step [log 043]
        if (ckpt_on) {
            const auto ckpt_now = std::chrono::steady_clock::now();
            const bool every_now = config.checkpoint_every_steps > 0 &&
                (step % static_cast<std::size_t>(config.checkpoint_every_steps) == 0);
            if (std::chrono::duration<double>(ckpt_now - ckpt_wall_last).count()
                    >= config.checkpoint_interval_s ||
                config.checkpoint_interval_s < 0 || every_now) {
                ckpt_write(step);
                ckpt_wall_last = std::chrono::steady_clock::now();
            }
        }

        if (step % 5000 == 0) {
            std::this_thread::yield();
        }

        // ===== periodic profile snapshot (instrumentation) =====
        if (rcp_profile_prints() && step % 2000 == 0) {
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
        if (rcp_profile_prints() && swap_on && swap_engine.n_sweeps_ > 0) {
            std::cerr << "  swap: sweeps=" << swap_engine.n_sweeps_
                      << " build=" << swap_engine.t_build_ << "s cand=" << swap_engine.t_cand_
                      << "s edges/sweep=" << (double)swap_engine.n_edges_total_ / swap_engine.n_sweeps_
                      << " cand/sweep=" << (double)swap_engine.n_cand_total_ / swap_engine.n_sweeps_
                      << " eval/sweep=" << (double)swap_engine.n_eval_total_ / swap_engine.n_sweeps_
                      << " accepted=" << swap_total_accepted
                      << " relayouts=" << swap_relayouts << "\n";
        }
        rcp_cycle_probe_dump();
        rcp_shadow_dump();
        rcp_sentinel_dump();
    }

    PackingResult result;
    // Cycle 14: convert flat x back to nested for the result API.
    // checkpoint v3: final write on loop exit (capped runs leave a resumable
    // blob). MUST precede ALL result assembly - several state arrays (D, box,
    // particle_origin, swap_log_buf) are std::move'd into the result below.
    if (ckpt_on && steps > ckpt_last_step_written) ckpt_write(steps);

    result.positions.assign(N, std::vector<double>(Ndim, 0.0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t d = 0; d < Ndim; ++d)
            result.positions[i][d] = x[i * Ndim + d];
    if (ell_on) {   // ellipsoids: back to the user's D (first axis); orientations out
        for (std::size_t i = 0; i < N; ++i) D[i] /= ell_bfac[i];
        if (ell_shape.kind != 0) {   // rods: length and radius per particle, in the output order
            result.rod_length.assign(N, 0.0); result.rod_radius.assign(N, 0.0);
            for (std::size_t i = 0; i < N; ++i) { result.rod_length[i] = ell_rodu[2 * i]; result.rod_radius[i] = ell_rodu[2 * i + 1]; }
        }
        result.orientations = ell_nested_orientations(ell_rot, N, Ndim);
        if (!input.aspect_ratios.empty()) {   // a table came in: the table goes out, in the output particle order
            result.aspect_ratios.assign(N, std::vector<double>(Ndim - 1, 1.0));
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t k = 1; k < Ndim; ++k) result.aspect_ratios[i][k - 1] = ell_q[i * Ndim + k];
        }
        result.torque_magnitude = ell_torque_magnitude(ell_torque, N, Ndim, Fmean, Lc);
        result.ell_permute_candidates = ell_perm_candidates;
        result.ell_capture_step = ell_cap_step; result.ell_capture_energy = ell_cap_energy;
        result.ell_capture_positions = std::move(ell_cap_x); result.ell_capture_diameters = std::move(ell_cap_D);
        result.ell_capture_orientations = std::move(ell_cap_rot); result.ell_capture_ratios = std::move(ell_cap_q);
        result.ell_capture_forces = std::move(ell_cap_F); result.ell_capture_torques = std::move(ell_cap_T);
        result.ell_capture_exchanges = std::move(ell_cap_exch);
        result.ell_capture_z = std::move(ell_cap_z);
    }
    result.diameters = std::move(D);
    result.particle_origin = std::move(particle_origin);
    result.box = std::move(box);
    result.walls = std::move(walls);
    // Trim the per-step history series to the number of steps actually
    // recorded. They are over-allocated to hist_size = step_cap + 1 for indexing
    // headroom, but the loop fills only [0, hist_filled); returning the full
    // allocation left a trailing zero, making every *_history read one element
    // longer than the step count. Return exactly the filled prefix.
    phi_history.resize(hist_filled);
    F_history.resize(hist_filled);
    U_history.resize(hist_filled);
    max_overlap_history.resize(hist_filled);
    mu_hist.resize(hist_filled);
    alpha_hist.resize(hist_filled);
    mu_flag_hist.resize(hist_filled);
    result.phi_history = std::move(phi_history);
    result.force_history = std::move(F_history);
    // Cycle 16: surface the diagnostic series to Python.
    result.max_overlap_history = std::move(max_overlap_history);
    result.mu_history = std::move(mu_hist);
    result.alpha_history = std::move(alpha_hist);
    result.mu_flag_history = std::move(mu_flag_hist);
    if (config.record_visibility) {
        zpos_hist.resize(hist_filled); mueff_hist.resize(hist_filled);
        result.zpos_history = std::move(zpos_hist);
        result.mu_eff_history = std::move(mueff_hist);
    }
    result.mu_hold_applied = mu_hold_applied;
    result.mu_ctrl_final_mult = mu_ctrl_mult;
    result.morph_log = std::move(morph_log_buf);       // [219]
    result.morph_lambda = morph_lambda;
    result.morph_updates = morph_updates;
    result.morph_holds = morph_holds;
    result.morph_stalled = morph_stalled ? 1 : 0;
    result.morph_rate_ceiling = morph_rate_ceiling;
    result.quench_residual = q_res; result.quench_global_residual = q_gres; result.quench_iters = q_iter;
    result.quench_converged = q_converged ? 1 : 0; result.quench_log = std::move(q_log); result.quench_held_force = q_held_force;
    if (quench_on && q_status == 0) q_status = 2;   // [225] the step loop ended before any stopping rule: report it as the cap
    result.quench_status = q_status;
    result.swap_geo_trades = static_cast<long>(swap_engine.geo_trades_);
    result.swap_geo_passes = static_cast<long>(swap_engine.geo_passes_);
    result.swap_geo_reopens = geo_reopens_done;
    result.swap_geo_cycle_trades = static_cast<long>(swap_engine.geo_cycle_trades_);
    result.swap_geo_cycle4_trades = static_cast<long>(swap_engine.geo_cycle4_trades_);   // [194]
    result.swap_geo_capk_trades = static_cast<long>(swap_engine.geo_capk_trades_);
    result.swap_geo_cage_trades = static_cast<long>(swap_engine.geo_cage_trades_);
    result.swap_geo_bottleneck_trades = static_cast<long>(swap_engine.geo_bottleneck_trades_);
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
    // SWAP diagnostics: the never-written reorder_log slot carries the per-sweep
    // event log (8 doubles per sweep) and, when the spike log is otherwise unused,
    // spike_log carries per-particle accepted-swap counts (N doubles). Both empty
    // when swap == 0, exactly as before.
    if (swap_on) {
        reorder_log_buf = std::move(swap_log_buf);
        if (!spike_log_enabled) {
            spike_log_buf.assign(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) spike_log_buf[k] = static_cast<double>(swap_engine.swap_count(k));
        }
    }
    result.spike_log = std::move(spike_log_buf);
    result.reorder_log = std::move(reorder_log_buf);
    result.swap_relayouts = swap_relayouts;
    result.swap_rescue_candidates = static_cast<std::size_t>(swap_engine.rescue_cand_total_);
    result.swap_rescue_accepted = static_cast<std::size_t>(swap_engine.rescue_acc_total_);
    result.swap_rescue_false_rej = static_cast<std::size_t>(swap_engine.rescue_frej_total_);
    result.swap_rescue_false_acc = static_cast<std::size_t>(swap_engine.rescue_facc_total_);
    result.swap_cycle3_candidates = static_cast<std::size_t>(swap_engine.cycle3_cand_total_);
    result.swap_cycle3_accepted = static_cast<std::size_t>(swap_engine.cycle3_acc_total_);
    result.swap_cycle3_blocked = static_cast<std::size_t>(swap_engine.cycle3_blocked_total_);
    result.swap_events = swap_engine.event_buf_;
    result.swap_event_width = static_cast<int>(swap_engine.event_width());
    (void)swap_total_accepted;
    // Report the number of COMPLETED packing iterations, which equals the
    // filled length of every per-step history series. The loop-entry counter
    // `steps` also counts a terminating check-only iteration (convergence hits
    // OscPhase::DONE and breaks before that iteration writes history), so it can
    // read one higher than the histories. Using hist_filled keeps the public
    // contract internally consistent: result.steps == len(phi/force/energy_history).
    result.resumed_from_step = ckpt_resumed_from;
    result.checkpoint_last_step = ckpt_last_step_written;
    result.steps = hist_filled;
    result.capped = (config.max_steps > 0 &&
                     hist_filled >= static_cast<std::size_t>(config.max_steps));
    result.phi = phi;
    result.max_min_dist = max_min_dist;
    result.force_magnitude = F_magnitude;

    const bool final_already_captured =
        capture_positions && steps > 0 && traj_last_captured == static_cast<std::size_t>(steps);
    const bool final_already_reported =
        observer && progress_interval > 0 && steps > 0 && (steps % progress_interval == 0);
    if ((capture_positions && !final_already_captured) || (observer && !final_already_reported)) {
        trace.steps.push_back(steps);
        trace.diameters.push_back(result.diameters);
        if (ell_on) trace.orientations.push_back(result.orientations);
        trace.phi.push_back(phi);
        trace.force.push_back(F_magnitude);
        trace.energy.push_back(U);
        trace.max_min_dist.push_back(max_min_dist);
        if (capture_positions) {
            trace.positions.push_back(result.positions);
            trace.origins.push_back(result.particle_origin);   // particle_origin itself was moved into the result above
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
