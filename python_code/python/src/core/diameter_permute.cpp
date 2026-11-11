#include "rcpgenerator/diameter_permute.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <chrono>

#include <omp.h>

namespace rcpgenerator {

namespace {
inline double e_term(double overlap) {
    return overlap > 0.0 ? 0.5 * overlap * overlap : 0.0;
}
}  // namespace

SwapEngine::SwapEngine(const SwapParams& p, std::uint64_t seed) : p_(p) {
    // splitmix64 state; the constant keeps seed 0 from producing a zero stream.
    rng_ = seed * 0x9E3779B97F4A7C15ULL + 0x53574150ULL;  // "SWAP"
    if (rng_ == 0) rng_ = 0x53574150ULL;
}

std::uint64_t SwapEngine::next_u64() {
    // splitmix64
    std::uint64_t z = (rng_ += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

double SwapEngine::next_unit() {
    return static_cast<double>(next_u64() >> 11) * (1.0 / 9007199254740992.0);
}

void SwapEngine::rebuild_ranks(const std::vector<double>& D0) {
    const std::size_t N = D0.size();
    at_rank_.resize(N);
    rank_of_.resize(N);
    std::iota(at_rank_.begin(), at_rank_.end(), 0u);
    std::stable_sort(at_rank_.begin(), at_rank_.end(),
                     [&D0](std::uint32_t a, std::uint32_t b) { return D0[a] < D0[b]; });
    for (std::size_t r = 0; r < N; ++r) rank_of_[at_rank_[r]] = static_cast<std::uint32_t>(r);
    rank_D0_.resize(N);
    for (std::size_t r = 0; r < N; ++r) rank_D0_[r] = D0[at_rank_[r]];
    if (p_.mode == 5) {
        if (p_.ctrl) {
            if (w_cur_ < 0.0) w_cur_ = p_.w_max;      // start wide; controller pulls down
            rebuild_window_table(w_cur_);
        } else if (p_.window_rel > 0.0) {
            rebuild_window_table(p_.window_rel);
        }
    }
    ranks_valid_ = true;
}

void SwapEngine::build_window_table_for(double w, std::vector<std::uint32_t>& wlo,
                                        std::vector<std::uint32_t>& whi) const {
    // for rank r, partners are ranks whose D0 lies within [D0(r)/(1+w), D0(r)*(1+w)]
    // (two-pointer sweep over the sorted rank_D0_, O(N))
    const std::size_t N = rank_D0_.size();
    wlo.resize(N); whi.resize(N);
    std::size_t lo = 0, hi = 0;
    for (std::size_t r = 0; r < N; ++r) {
        const double d = rank_D0_[r];
        while (lo < r && rank_D0_[lo] * (1.0 + w) < d) ++lo;
        if (hi < r) hi = r;
        while (hi + 1 < N && rank_D0_[hi + 1] <= d * (1.0 + w)) ++hi;
        wlo[r] = static_cast<std::uint32_t>(lo);
        whi[r] = static_cast<std::uint32_t>(hi);
    }
}

void SwapEngine::rebuild_window_table(double w) {
    build_window_table_for(w, win_lo_, win_hi_);
    wmax_cached_ = 0;   // banded sweep re-derives the max window span lazily
    rwin_built_ = false;  // sorted D0 may have been rebuilt with the ranks
}

// ---- checkpoint support [log 043] -----------------------------------------------
namespace {
constexpr std::uint32_t kSwapStateVersion = 2;
bool wr(std::FILE* fp, const void* p, std::size_t n) {
    return std::fwrite(p, 1, n, fp) == n;
}
bool rd(std::FILE* fp, void* p, std::size_t n) {
    return std::fread(p, 1, n, fp) == n;
}
template <class T>
bool wr_vec(std::FILE* fp, const std::vector<T>& v) {
    const std::uint64_t n = v.size();
    if (!wr(fp, &n, 8)) return false;
    return n == 0 || wr(fp, v.data(), n * sizeof(T));
}
template <class T>
bool rd_vec(std::FILE* fp, std::vector<T>& v) {
    std::uint64_t n = 0;
    if (!rd(fp, &n, 8)) return false;
    v.resize(static_cast<std::size_t>(n));
    return n == 0 || rd(fp, v.data(), static_cast<std::size_t>(n) * sizeof(T));
}
}  // namespace

bool SwapEngine::save_state(std::FILE* fp) const {
    bool ok = wr(fp, &kSwapStateVersion, 4);
    const std::uint8_t rv = ranks_valid_ ? 1 : 0, rwb = rwin_built_ ? 1 : 0;
    const std::uint64_t wmax64 = wmax_cached_;
    ok = ok && wr(fp, &rng_, 8) && wr(fp, &rv, 1) && wr(fp, &rwb, 1)
            && wr(fp, &w_cur_, 8) && wr(fp, &wmax64, 8)
            && wr(fp, &t_build_, 8) && wr(fp, &t_cand_, 8)
            && wr(fp, &n_edges_total_, 8) && wr(fp, &n_eval_total_, 8)
            && wr(fp, &n_cand_total_, 8) && wr(fp, &n_sweeps_, 8)
            && wr(fp, &rescue_cand_total_, 8) && wr(fp, &rescue_acc_total_, 8)
            && wr(fp, &rescue_frej_total_, 8) && wr(fp, &rescue_facc_total_, 8)
            && wr(fp, &cycle3_cand_total_, 8) && wr(fp, &cycle3_acc_total_, 8)
            && wr(fp, &cycle3_blocked_total_, 8);
    ok = ok && wr_vec(fp, event_buf_);
    ok = ok && wr_vec(fp, rank_of_) && wr_vec(fp, at_rank_)
            && wr_vec(fp, rank_D0_) && wr_vec(fp, win_lo_) && wr_vec(fp, win_hi_)
            && wr_vec(fp, rwin_lo_) && wr_vec(fp, rwin_hi_)
            && wr_vec(fp, swap_count_);
    return ok;
}

bool SwapEngine::load_state(std::FILE* fp) {
    std::uint32_t ver = 0;
    if (!rd(fp, &ver, 4) || ver != kSwapStateVersion) return false;
    std::uint8_t rv = 0, rwb = 0;
    std::uint64_t wmax64 = 0;
    bool ok = rd(fp, &rng_, 8) && rd(fp, &rv, 1) && rd(fp, &rwb, 1)
           && rd(fp, &w_cur_, 8) && rd(fp, &wmax64, 8)
           && rd(fp, &t_build_, 8) && rd(fp, &t_cand_, 8)
           && rd(fp, &n_edges_total_, 8) && rd(fp, &n_eval_total_, 8)
           && rd(fp, &n_cand_total_, 8) && rd(fp, &n_sweeps_, 8)
           && rd(fp, &rescue_cand_total_, 8) && rd(fp, &rescue_acc_total_, 8)
           && rd(fp, &rescue_frej_total_, 8) && rd(fp, &rescue_facc_total_, 8)
           && rd(fp, &cycle3_cand_total_, 8) && rd(fp, &cycle3_acc_total_, 8)
           && rd(fp, &cycle3_blocked_total_, 8);
    ok = ok && rd_vec(fp, event_buf_);
    ok = ok && rd_vec(fp, rank_of_) && rd_vec(fp, at_rank_)
            && rd_vec(fp, rank_D0_) && rd_vec(fp, win_lo_) && rd_vec(fp, win_hi_)
            && rd_vec(fp, rwin_lo_) && rd_vec(fp, rwin_hi_)
            && rd_vec(fp, swap_count_);
    if (ok) {
        ranks_valid_ = (rv != 0);
        rwin_built_ = (rwb != 0);
        wmax_cached_ = static_cast<std::size_t>(wmax64);
    }
    return ok;
}

// Surgical injection (ME0-ii [log 054]): forced trade with sweep-identical
// bookkeeping. relief is logged as NaN (no evaluation happens); the twin
// protocol measures the causal effect instead.
bool SwapEngine::inject(std::size_t N, std::vector<double>& D,
                        std::vector<double>& D0,
                        std::vector<std::uint32_t>& touched,
                        std::uint32_t i, std::uint32_t j)
{
    if (i == j || i >= N || j >= N) return false;
    if (!ranks_valid_) rebuild_ranks(D0);
    const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
    const double d_rel = (D0[i] > D0[j] ? D0[i] / D0[j] : D0[j] / D0[i]) - 1.0;
    const double d_abs = 0.5 * std::fabs(D[i] - D[j]);
    std::swap(D[i], D[j]);
    std::swap(D0[i], D0[j]);
    rank_of_[i] = qj; rank_of_[j] = qi;
    at_rank_[qi] = j; at_rank_[qj] = i;
    touched.push_back(i);
    touched.push_back(j);
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    ++swap_count_[i]; ++swap_count_[j];
    const std::size_t EW = event_width();
    if (p_.event_log_stride > 0 && !event_log_full_ &&
        event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
        event_buf_.push_back(static_cast<double>(current_step_));
        event_buf_.push_back(static_cast<double>(i));
        event_buf_.push_back(static_cast<double>(j));
        event_buf_.push_back(d_rel);
        event_buf_.push_back(d_abs);
        event_buf_.push_back(std::numeric_limits<double>::quiet_NaN());
        event_buf_.push_back(-1.0);
        event_buf_.push_back(-1.0);
        for (std::size_t q = 8; q < EW; ++q)
            event_buf_.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    return true;
}

// ME0-iv [log 067]: per-particle cage radius and active-contact fraction from
// the near-pair edges of this sweep step (positions and radii as the forces saw
// them). Only ever called with event_log_geom on; no effect on the trade path.
void SwapEngine::geom_from_edges(std::size_t N, const std::vector<double>& D,
                                 const SwapEdgeSink& edges,
                                 const std::vector<std::size_t>& z)
{
    // Neighbour shell: r < (1 + kGeomTol) * (R_i + R_k), a FIXED relative
    // tolerance (the collection margin itself tracks the controller window and
    // would make the shell w-dependent). rcp_generator widens the collection
    // margin to cover this shell whenever geom is on.
    constexpr double kGeomTol = 0.10;
    std::vector<double> csum(N, 0.0);
    std::vector<std::uint32_t> nn(N, 0u);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            const double s = 0.5 * (D[ed.i] + D[ed.j]);
            if (r >= (1.0 + kGeomTol) * s) continue;
            csum[ed.i] += r - 0.5 * D[ed.j]; ++nn[ed.i];
            csum[ed.j] += r - 0.5 * D[ed.i]; ++nn[ed.j];
        }
    }
    cage_.assign(N, std::numeric_limits<double>::quiet_NaN());
    fact_.assign(N, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < N; ++i) {
        if (nn[i]) {
            cage_[i] = csum[i] / static_cast<double>(nn[i]);
            const double zi = i < z.size() ? static_cast<double>(z[i]) : 0.0;
            fact_[i] = std::min(1.0, zi / static_cast<double>(nn[i]));
        }
    }
}

// Stage 3 [log 082]: one geometric sort pass (see SwapParams::geo). Positions
// are fixed; the inscribed radius c_i = min_k (r_ik - R_k) over the 1.5x shell
// comes from this step's near-pair edges (rcp_generator widens the collection
// margin to 0.5 D_max on geo steps). No RNG is consumed at T = 0.
std::uint64_t SwapEngine::geo_pass(std::size_t N, std::vector<double>& D,
                                   std::vector<double>& D0,
                                   const SwapEdgeSink& edges,
                                   std::vector<std::uint32_t>& touched)
{
    if (N < 2 || !edges.want_edges) return 0;
    if (!ranks_valid_) rebuild_ranks(D0);
    // [135] annealed temperature: fixed geo_T unless geo_T_decay<1 & geo_T_end>=0.
    const bool anneal = (p_.geo_T_decay < 1.0) && (p_.geo_T_end >= 0.0);
    if (anneal && geo_T_cur_ < 0.0) geo_T_cur_ = p_.geo_T;          // lazy init on first pass
    const double Teff = anneal ? geo_T_cur_ : p_.geo_T;
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> cmin(N, INF);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            cmin[ed.i] = std::min(cmin[ed.i], r - 0.5 * D[ed.j]);
            cmin[ed.j] = std::min(cmin[ed.j], r - 0.5 * D[ed.i]);
        }
    }
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    struct Cand { double gain; std::uint32_t i, j; };
    std::vector<Cand> cands;
    cands.reserve(N * 8);
    const std::uint32_t* at_rank = at_rank_.data();
    const int cap_p = std::max(1, p_.geo_max_partners);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!(cmin[i] < INF)) continue;
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0) continue;
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap_p - 1) / cap_p);
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!(cmin[j] < INF)) continue;
            const double dR = 0.5 * (D[j] - D[i]);
            if (dR <= 0.0) continue;
            const double gain = 2.0 * dR * (cmin[i] - cmin[j]);   // = -dM; > 0 improves
            if (Teff <= 0.0 && gain <= 0.0) continue;
            cands.push_back({gain, i, j});
        }
    }
    // ---- H19 [172]: non-separable rescore of the top candidates -------------------
    // Replace each top candidate's separable gain with the TRUE dJ' of the concave-capped
    // objective J' = sum_p min(cmin_p, ct), summed over N(i) U N(j) (the third-party cage
    // shifts the geo gain drops). Frozen membership + frozen positions = C5-clean.
    const bool nonsep_in_window =                                       // [190] H26 window
        static_cast<long>(current_step_) >= p_.geo_nonsep_start &&
        (p_.geo_nonsep_stop < 0 || static_cast<long>(current_step_) < p_.geo_nonsep_stop);
    const bool nonsep_on = (p_.geo_nonsep == 1 && nonsep_in_window) ||
                           (p_.geo_nonsep == 2 && geo_reopen_active_);   // 2 = reopen-only
    if (nonsep_on && !cands.empty()) {
        std::sort(cands.begin(), cands.end(),
                  [](const Cand& a, const Cand& b) { return a.gain > b.gain; });
        const std::size_t topK = std::min(cands.size(),
            static_cast<std::size_t>(std::max(1L, static_cast<long>(p_.geo_nsk))
                                     * std::max(1L, p_.geo_B)));
        cands.resize(topK);
        // CSR adjacency from the near-pair edges (same pass-start cutoff as cmin).
        std::vector<std::uint32_t> deg(N, 0);
        for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
            const SwapEdge* v = edges.per_thread[t].data();
            const std::size_t n = edges.count[t];
            for (std::size_t q = 0; q < n; ++q) {
                const SwapEdge& ed = v[q];
                const double r = std::sqrt(ed.d2);
                if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
                ++deg[ed.i]; ++deg[ed.j];
            }
        }
        std::vector<std::size_t> off(N + 1, 0);
        for (std::size_t p = 0; p < N; ++p) off[p + 1] = off[p] + deg[p];
        const std::size_t NE = off[N];
        std::vector<std::uint32_t> adj(NE);
        std::vector<double> adr(NE);
        std::vector<std::size_t> cur(off.begin(), off.end() - 1);
        for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
            const SwapEdge* v = edges.per_thread[t].data();
            const std::size_t n = edges.count[t];
            for (std::size_t q = 0; q < n; ++q) {
                const SwapEdge& ed = v[q];
                const double r = std::sqrt(ed.d2);
                if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
                adj[cur[ed.i]] = ed.j; adr[cur[ed.i]] = r; ++cur[ed.i];
                adj[cur[ed.j]] = ed.i; adr[cur[ed.j]] = r; ++cur[ed.j];
            }
        }
        // Non-separable objective J' = sum_p D_p * g(cmin_p), the D-WEIGHTED (big-D wants
        // room, like separable geo) but CONCAVE (g=sqrt: tight cages weigh more) and LIVE-
        // cmin (coupled to neighbours' diameters) form. geo_ctarget is a small additive
        // clearance offset c0 (units of mean R) so g is smooth through 0; concavity (hence
        // the non-separability that beats 2-opt) is strongest for the tightest cages.
        double meanD = 0.0; for (std::size_t p = 0; p < N; ++p) meanD += D[p];
        meanD /= static_cast<double>(N);
        const double c0 = p_.geo_ctarget * 0.5 * meanD;
        auto g = [c0](double c) { const double a = c + c0; return a > 0.0 ? std::sqrt(a) : 0.0; };
        // cmin at site p with the diameters at sites si,sj swapped (membership frozen).
        auto cmin_sub = [&](std::uint32_t p, std::uint32_t si, std::uint32_t sj,
                            double Dsi, double Dsj) -> double {
            double m = INF;
            for (std::size_t e = off[p]; e < off[p + 1]; ++e) {
                const std::uint32_t k = adj[e];
                double Dk = D[k];
                if (k == si) Dk = Dsi; else if (k == sj) Dk = Dsj;
                const double val = adr[e] - 0.5 * Dk;
                if (val < m) m = val;
            }
            return m;
        };
        std::vector<std::uint32_t> stamp(N, 0);
        std::uint32_t tok = 0;
        std::vector<std::uint32_t> aff; aff.reserve(64);
        for (Cand& c : cands) {
            const std::uint32_t i = c.i, j = c.j;
            const double Di = D[i], Dj = D[j];
            // affected = {i,j} (their D flips) U N(i) U N(j) (their cmin shifts)
            ++tok; aff.clear();
            stamp[i] = tok; aff.push_back(i);
            stamp[j] = tok; aff.push_back(j);
            for (std::size_t e = off[i]; e < off[i + 1]; ++e) {
                const std::uint32_t p = adj[e];
                if (stamp[p] != tok) { stamp[p] = tok; aff.push_back(p); }
            }
            for (std::size_t e = off[j]; e < off[j + 1]; ++e) {
                const std::uint32_t p = adj[e];
                if (stamp[p] != tok) { stamp[p] = tok; aff.push_back(p); }
            }
            double dJ = 0.0;
            for (const std::uint32_t p : aff) {
                const double c_old = (cmin[p] < INF) ? cmin[p] : 0.0;
                const double cs = cmin_sub(p, i, j, Dj, Di);   // site i gets Dj, j gets Di
                const double c_new = (cs < INF) ? cs : 0.0;
                const double D_new = (p == i) ? Dj : (p == j) ? Di : D[p];
                dJ += D_new * g(c_new) - D[p] * g(c_old);
            }
            c.gain = dJ;
        }
        if (Teff <= 0.0) {
            cands.erase(std::remove_if(cands.begin(), cands.end(),
                        [](const Cand& c) { return !(c.gain > 0.0); }), cands.end());
        }
    }
    if (Teff <= 0.0) {
        std::sort(cands.begin(), cands.end(),
                  [](const Cand& a, const Cand& b) { return a.gain > b.gain; });
    } else {
        for (std::size_t k = cands.size(); k > 1; --k) {      // Fisher-Yates
            const std::size_t m = static_cast<std::size_t>(next_u64() % k);
            std::swap(cands[k - 1], cands[m]);
        }
    }
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Cand& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j]) continue;
        if (Teff > 0.0 && c.gain <= 0.0 &&
            next_unit() >= std::exp(c.gain / Teff)) continue;
        used[c.i] = used[c.j] = 1;
        const std::uint32_t i = c.i, j = c.j;          // i grows (smaller D0)
        const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
        const double d_rel = D0[j] / D0[i] - 1.0;
        const double dR = 0.5 * (D[j] - D[i]);
        std::swap(D[i], D[j]);
        std::swap(D0[i], D0[j]);
        rank_of_[i] = qj; rank_of_[j] = qi;
        at_rank_[qi] = j; at_rank_[qj] = i;
        touched.push_back(i); touched.push_back(j);
        ++swap_count_[i]; ++swap_count_[j];
        ++acc;
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back(d_rel);
            event_buf_.push_back(dR);
            event_buf_.push_back(c.gain);          // misfit gain (-dM), not an energy relief
            event_buf_.push_back(-2.0);            // marker: geometric pass
            event_buf_.push_back(-2.0);
            for (std::size_t q = 8; q < EW; ++q)
                event_buf_.push_back(q == 8 ? cmin[i] : q == 9 ? cmin[j] : nan);
        }
    }
    geo_trades_ += acc; ++geo_passes_;
    if (anneal) geo_T_cur_ = std::max(p_.geo_T_end, geo_T_cur_ * p_.geo_T_decay);  // [135] cool
    return acc;
}

// H7 [121/122]: cage-capacity sort pass. Runs AFTER geo_pass on the same geo step.
// c_min = cap[0] (single facet) is degenerate at a jammed terminal; the cage that
// actually constrains a sphere in 3D is set by >= 4 neighbours. We keep the K
// smallest per-neighbour capacities cap[i][0..K-1] (K = geo_capk) and re-rank D
// against a cage signal (geo_capmode): 1 = K-th smallest; 2 = mean of K; 3 = tiered
// "A else B" - take the c_min(cap[0]) gain if it improves, ELSE the cage-composite
// (mean-of-K) gain (this taps the residual that pairwise-c_min sorting is blind to).
// OVERLAP-TOLERANT by design: no feasibility reject; transient overlaps are resolved
// by the solver's relaxation into the cage room (user reframe [121]). T=0, exact
// descent, each particle once, <= geo_B swaps. Objective J = sum_p D_p * cap_p.
std::uint64_t SwapEngine::geo_capk_pass(std::size_t N, std::vector<double>& D,
                                        std::vector<double>& D0,
                                        const SwapEdgeSink& edges,
                                        std::vector<std::uint32_t>& touched)
{
    if (N < 2 || !edges.want_edges) return 0;
    if (!ranks_valid_) rebuild_ranks(D0);
    const double INF = std::numeric_limits<double>::infinity();
    const int K = std::min(8, std::max(1, p_.geo_capk));
    const int mode = p_.geo_capmode;
    std::vector<double> cap(N * static_cast<std::size_t>(K), INF);   // K smallest, ascending
    auto insert_cap = [&](std::uint32_t p, double val) {
        double* c = &cap[static_cast<std::size_t>(p) * K];
        if (!(val < c[K - 1])) return;                 // not among the K smallest
        int m = K - 1;
        while (m > 0 && val < c[m - 1]) { c[m] = c[m - 1]; --m; }
        c[m] = val;
    };
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            insert_cap(ed.i, r - 0.5 * D[ed.j]);
            insert_cap(ed.j, r - 0.5 * D[ed.i]);
        }
    }
    // per-particle signals: cmin = cap[0]; composite = mean of the K smallest (finite).
    std::vector<double> cmin(N), comp(N);
    std::vector<unsigned char> ok(N, 0);
    for (std::size_t p = 0; p < N; ++p) {
        const double* c = &cap[p * K];
        if (!(c[0] < INF)) continue;                   // no neighbours at all
        cmin[p] = c[0];
        double sum = 0.0; int cnt = 0;
        for (int m = 0; m < K; ++m) if (c[m] < INF) { sum += c[m]; ++cnt; }
        // mode 1 uses the K-th smallest but only when a full cage (K finite) exists
        comp[p] = (mode == 1) ? (cnt == K ? c[K - 1] : INF) : (sum / cnt);
        ok[p] = 1;
    }
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    struct Cand { double gain; std::uint32_t i, j; };
    std::vector<Cand> cands;
    cands.reserve(N * 4);
    const std::uint32_t* at_rank = at_rank_.data();
    const int cap_p = std::max(1, p_.geo_max_partners);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!ok[i]) continue;
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0) continue;
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap_p - 1) / cap_p);
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!ok[j]) continue;
            const double dR = 0.5 * (D[j] - D[i]);
            if (dR <= 0.0) continue;                    // i grows (smaller D0)
            double gain;
            if (mode == 3) {                            // tiered: c_min first, else cage
                const double gA = 2.0 * dR * (cmin[i] - cmin[j]);
                if (gA > 0.0) gain = gA;
                else {
                    if (!(comp[i] < INF) || !(comp[j] < INF)) continue;
                    gain = 2.0 * dR * (comp[i] - comp[j]);   // the residual channel
                }
            } else {                                    // 1 = K-th, 2 = mean-of-K
                if (!(comp[i] < INF) || !(comp[j] < INF)) continue;
                gain = 2.0 * dR * (comp[i] - comp[j]);
            }
            if (gain <= 0.0) continue;
            cands.push_back({gain, i, j});
        }
    }
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) { return a.gain > b.gain; });
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Cand& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j]) continue;
        used[c.i] = used[c.j] = 1;
        const std::uint32_t i = c.i, j = c.j;
        const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
        const double d_rel = D0[j] / D0[i] - 1.0;
        const double dR = 0.5 * (D[j] - D[i]);
        std::swap(D[i], D[j]);
        std::swap(D0[i], D0[j]);
        rank_of_[i] = qj; rank_of_[j] = qi;
        at_rank_[qi] = j; at_rank_[qj] = i;
        touched.push_back(i); touched.push_back(j);
        ++swap_count_[i]; ++swap_count_[j];
        ++acc;
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back(d_rel);
            event_buf_.push_back(dR);
            event_buf_.push_back(c.gain);
            event_buf_.push_back(-4.0);            // marker: cage-capacity pass
            event_buf_.push_back(-4.0);
            for (std::size_t q = 8; q < EW; ++q)
                event_buf_.push_back(q == 8 ? cmin[i] : q == 9 ? cmin[j] : nan);
        }
    }
    geo_capk_trades_ += acc;
    return acc;
}

// H8 [129/130/131]: true-cage OVERLAP-ENERGY pairwise pass. Runs AFTER geo_pass on
// the same geo step. The move runs above jamming under compression; the overlap is
// the relaxation signal (user). Overlap of a trial diameter D at site i with contact
// k is max(0, D/2 - gap_ik), gap_ik = r_ik - R_k (depends on the NEIGHBOUR radius,
// not on i's diameter => frozen at pass start like geo's c_min). Cage overlap energy
// E(D@i) = sum_k max(0, D/2 - gap_ik)^2 over i's CAGE_K nearest contacts (>= any 3D
// contact count). A rank-window pair (i,j) is swapped iff it lowers total energy:
// E(Dj@i)+E(Di@j) < E(Di@i)+E(Dj@j). Greedy (T=0), each particle once, <= geo_B
// swaps. OVERLAP-TOLERANT (this IS the overlap objective; no feasibility reject).
// geo_cage_trades_ over time is the pairwise-completeness certificate: -> 0 means no
// overlap-improving pairwise swap remains (2-opt-optimal on the true objective) -
// which c_min (a linear proxy) cannot certify. Scored at the CURRENT position; the
// static probe [130] showed recentering to the tetra centre adds <2.5% (and the user
// does not want a position move here). No RNG consumed.
// [133] DIAGNOSTIC deepest-point (cage centre) for the position-move variant. Ports
// the smoke-tested analysis/cage_true.py solver: BACKTRACKING Chebyshev ascent of
// min_k(|p-x_k|-R_k) - move away from ALL near-binding walls at once (normalized sum
// of outward unit vectors over the active set; = 0 exactly when they surround p =>
// converged), accept a step only if it raises the min-clearance (so p cannot jump a
// saddle into a neighbouring cell), displacement-capped at max_disp. Returns the
// centre into out[0..Ndim-1]. m = number of cage neighbours provided.
// Returns the SURROUND quality at the found centre: resultant of the active outward
// unit vectors / n_active, in [0,1] (0 = fully surrounds, 1 = one-sided). The caller
// must move ONLY when this is small (< 0.5) - else the "centre" is on a void wall and
// moving there pushes the particle into unlisted neighbours (destructive).
static double cage_center_dp(const double* xi, const double* xk, const double* Rk,
                             int m, const double* box, int Ndim, double max_disp,
                             double* out)
{
    auto mimg = [&](double d, int c) { return d - box[c] * std::round(d / box[c]); };
    double Rmed;
    {   // median of Rk (small m; simple copy-sort)
        double tmp[32]; int mm = m < 32 ? m : 32;
        for (int k = 0; k < mm; ++k) tmp[k] = Rk[k];
        std::sort(tmp, tmp + mm); Rmed = tmp[mm / 2];
    }
    double p[3]; for (int c = 0; c < Ndim; ++c) p[c] = xi[c];
    auto minclear = [&](const double* q) {
        double best = std::numeric_limits<double>::infinity();
        for (int k = 0; k < m; ++k) {
            double d2 = 0.0;
            for (int c = 0; c < Ndim; ++c) { double dd = mimg(q[c] - xk[k * Ndim + c], c); d2 += dd * dd; }
            double g = std::sqrt(d2) - Rk[k];
            if (g < best) best = g;
        }
        return best;
    };
    double best = minclear(p);
    double step = 0.25 * Rmed;
    for (int it = 0; it < 160; ++it) {
        // active set = neighbours within 0.15*Rmed of the current min; ascent dir =
        // normalized sum of outward units (from each active wall toward p).
        double dir[3] = {0, 0, 0}; int nact = 0;
        for (int k = 0; k < m; ++k) {
            double dv[3], d2 = 0.0;
            for (int c = 0; c < Ndim; ++c) { dv[c] = mimg(p[c] - xk[k * Ndim + c], c); d2 += dv[c] * dv[c]; }
            double dist = std::sqrt(d2), g = dist - Rk[k];
            if (g < best + 0.15 * Rmed) {
                for (int c = 0; c < Ndim; ++c) dir[c] += dv[c] / (dist + 1e-12);
                ++nact;
            }
        }
        double nd = 0.0; for (int c = 0; c < Ndim; ++c) nd += dir[c] * dir[c]; nd = std::sqrt(nd);
        if (nd < 1e-9) break;                              // surrounded => at the centre
        double pn[3];
        for (int c = 0; c < Ndim; ++c) pn[c] = p[c] + step * dir[c] / nd;
        double disp = 0.0;
        for (int c = 0; c < Ndim; ++c) { double e = pn[c] - xi[c]; disp += e * e; }
        if (std::sqrt(disp) > max_disp) { step *= 0.5; if (step < 1e-4 * Rmed) break; continue; }
        double gn = minclear(pn);
        if (gn > best + 1e-12) { for (int c = 0; c < Ndim; ++c) p[c] = pn[c]; best = gn; }
        else { step *= 0.5; if (step < 1e-4 * Rmed) break; }
    }
    for (int c = 0; c < Ndim; ++c) out[c] = p[c];
    // surround quality at the final point: resultant of active outward units / n_active
    double res[3] = {0, 0, 0}; int nact = 0;
    for (int k = 0; k < m; ++k) {
        double dv[3], d2 = 0.0;
        for (int c = 0; c < Ndim; ++c) { dv[c] = mimg(p[c] - xk[k * Ndim + c], c); d2 += dv[c] * dv[c]; }
        double dist = std::sqrt(d2), g = dist - Rk[k];
        if (g < best + 0.05 * Rmed) { for (int c = 0; c < Ndim; ++c) res[c] += dv[c] / (dist + 1e-12); ++nact; }
    }
    if (nact < 4) return 1.0;                               // under-coordinated => not a cage
    double rn = 0.0; for (int c = 0; c < Ndim; ++c) rn += res[c] * res[c];
    return std::sqrt(rn) / nact;
}

std::uint64_t SwapEngine::geo_cage_pass(std::size_t N, std::vector<double>& D,
                                        std::vector<double>& D0,
                                        const SwapEdgeSink& edges,
                                        std::vector<std::uint32_t>& touched,
                                        double* x, const double* box, int Ndim)
{
    if (N < 2 || !edges.want_edges) return 0;
    const bool do_move = (p_.geo_cage_move != 0) && (x != nullptr) && (box != nullptr);
    if (!ranks_valid_) rebuild_ranks(D0);
    const double INF = std::numeric_limits<double>::infinity();
    const int K = 12;                              // cage contacts that can overlap (>= 3D max)
    std::vector<double> gap(N * static_cast<std::size_t>(K), INF);   // K smallest gaps, ascending
    // parallel neighbour indices (only when the diagnostic position move needs them)
    std::vector<std::uint32_t> gidx(do_move ? N * static_cast<std::size_t>(K) : 0, 0);
    auto insert_gap = [&](std::uint32_t p, double val, std::uint32_t nbr) {
        double* c = &gap[static_cast<std::size_t>(p) * K];
        if (!(val < c[K - 1])) return;
        int m = K - 1;
        if (do_move) {
            std::uint32_t* gi = &gidx[static_cast<std::size_t>(p) * K];
            while (m > 0 && val < c[m - 1]) { c[m] = c[m - 1]; gi[m] = gi[m - 1]; --m; }
            c[m] = val; gi[m] = nbr;
        } else {
            while (m > 0 && val < c[m - 1]) { c[m] = c[m - 1]; --m; }
            c[m] = val;
        }
    };
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            insert_gap(ed.i, r - 0.5 * D[ed.j], ed.j);   // gap_ik uses NEIGHBOUR radius R_k
            insert_gap(ed.j, r - 0.5 * D[ed.i], ed.i);
        }
    }
    std::vector<unsigned char> ok(N, 0);
    for (std::size_t p = 0; p < N; ++p) ok[p] = (gap[p * K] < INF) ? 1 : 0;
    // E(diam @ site p) = sum_k max(0, diam/2 - gap_pk)^2 over finite gaps.
    auto energy = [&](double diam, std::uint32_t p) {
        const double R = 0.5 * diam;
        const double* c = &gap[static_cast<std::size_t>(p) * K];
        double e = 0.0;
        for (int m = 0; m < K; ++m) {
            if (!(c[m] < INF)) break;
            const double ov = R - c[m];
            if (ov > 0.0) e += ov * ov;
        }
        return e;
    };
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    struct Cand { double gain; std::uint32_t i, j; };
    std::vector<Cand> cands;
    cands.reserve(N * 4);
    const std::uint32_t* at_rank = at_rank_.data();
    const int cap_p = std::max(1, p_.geo_max_partners);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!ok[i]) continue;
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0) continue;
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap_p - 1) / cap_p);
        const double Di = D[i];
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!ok[j]) continue;
            const double Dj = D[j];
            if (Di == Dj) continue;
            // gain = keep - swap; > 0 => swapping diameters lowers total overlap energy
            const double gain = (energy(Di, i) + energy(Dj, j))
                              - (energy(Dj, i) + energy(Di, j));
            if (gain <= 0.0) continue;
            cands.push_back({gain, i, j});
        }
    }
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) { return a.gain > b.gain; });
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Cand& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j]) continue;
        used[c.i] = used[c.j] = 1;
        const std::uint32_t i = c.i, j = c.j;
        const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
        const double d_rel = D0[j] / D0[i] - 1.0;
        const double dR = 0.5 * (D[j] - D[i]);
        std::swap(D[i], D[j]);
        std::swap(D0[i], D0[j]);
        rank_of_[i] = qj; rank_of_[j] = qi;
        at_rank_[qi] = j; at_rank_[qj] = i;
        touched.push_back(i); touched.push_back(j);
        ++swap_count_[i]; ++swap_count_[j];
        ++acc;
        if (do_move) {
            // DIAGNOSTIC [133]: place each swapped particle at its cage centre
            // (deepest point of its cell). Position move = OUTSIDE C5; phi from a
            // move-enabled run is an insight number, never a milestone. Gather the
            // (<=K) cage neighbours' positions/radii and solve the Chebyshev centre;
            // cap the displacement at the particle's post-swap radius.
            double xk[3 * 12], Rk[12], ctr[3];
            for (std::uint32_t pmv : {i, j}) {
                const std::uint32_t* gi = &gidx[static_cast<std::size_t>(pmv) * K];
                const double* gv = &gap[static_cast<std::size_t>(pmv) * K];
                int m = 0;
                for (int kk = 0; kk < K; ++kk) {
                    if (!(gv[kk] < INF)) break;
                    const std::uint32_t nb = gi[kk];
                    for (int cc = 0; cc < Ndim; ++cc) xk[m * Ndim + cc] = x[static_cast<std::size_t>(nb) * Ndim + cc];
                    Rk[m] = 0.5 * D[nb];
                    ++m;
                }
                if (m < 4) continue;                        // no enclosing cage -> don't move
                const double surround = cage_center_dp(&x[static_cast<std::size_t>(pmv) * Ndim],
                                                       xk, Rk, m, box, Ndim, 0.5 * D[pmv], ctr);
                if (surround >= 0.5) continue;              // one-sided => not a real cage; don't move
                for (int cc = 0; cc < Ndim; ++cc) {
                    double v = ctr[cc];
                    v -= box[cc] * std::floor(v / box[cc]);  // wrap into [0, box)
                    x[static_cast<std::size_t>(pmv) * Ndim + cc] = v;
                }
            }
        }
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back(d_rel);
            event_buf_.push_back(dR);
            event_buf_.push_back(c.gain);          // overlap-energy relief (keep - swap)
            event_buf_.push_back(do_move ? -6.0 : -5.0);  // -5 diam-only, -6 with position move
            event_buf_.push_back(do_move ? -6.0 : -5.0);
            for (std::size_t q = 8; q < EW; ++q)
                event_buf_.push_back(q == 8 ? gap[static_cast<std::size_t>(i) * K]
                                            : q == 9 ? gap[static_cast<std::size_t>(j) * K] : nan);
        }
    }
    geo_cage_trades_ += acc;
    return acc;
}

// H12 [146/148]: BOTTLENECK-RELIEF pass (user's "final squeeze"). Relieves the
// tightest contacts: for a bottleneck site i (small kappa_i = c_i/R_i), swap in a
// SMALLER diameter from a partner j whose site can ACCOMMODATE i's bigger diameter
// (c_j >= R_i, i.e. kappa_j stays >= 1 after taking R_i). Accept iff the swap raises
// the pair's min kappa (relieves the bottleneck without creating a worse one). c_i
// frozen at pass start (as geo_pass); the gentle small-mu relax between passes (the
// reopen loop) refreshes c and regenerates supply. Greedy T=0, each particle once,
// <= geo_B. Feasibility-RESPECTING (accommodate condition) - unlike the overlap-
// tolerant cage pass. Objective: max increase of the pair bottleneck min(kappa).
std::uint64_t SwapEngine::geo_bottleneck_pass(std::size_t N, std::vector<double>& D,
                                              std::vector<double>& D0,
                                              const SwapEdgeSink& edges,
                                              std::vector<std::uint32_t>& touched)
{
    if (N < 2 || !edges.want_edges) return 0;
    if (!ranks_valid_) rebuild_ranks(D0);
    const double INF = std::numeric_limits<double>::infinity();

    // -----------------------------------------------------------------------
    // H12c [156] mode 3: the user's EXACT algorithm — shrink to zero overlap
    // (one just-touching pair), relieve THAT pair by swapping its larger
    // diameter to a loose site that accommodates it, re-inflate (recompute the
    // just-touching scale s*), repeat until no swap raises s*. Positions are
    // fixed; s* is an evaluation device (the run's kappa is untouched). The
    // rung's small-mu ADAM steps between passes are the outer relax.
    if (p_.geo_bottleneck == 3) {
        // flat near-edge list (i, j, r) at the current geometry
        struct FE { std::uint32_t i, j; double r; };
        std::vector<FE> Eg;
        Eg.reserve(N * 6);
        for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
            const SwapEdge* v = edges.per_thread[t].data();
            const std::size_t n = edges.count[t];
            for (std::size_t q = 0; q < n; ++q) {
                const SwapEdge& ed = v[q];
                const double r = std::sqrt(ed.d2);
                if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
                Eg.push_back({ed.i, ed.j, r});
            }
        }
        if (Eg.empty()) return 0;
        // s* = min over edges of r / (Ri+Rj); tightest pair = the argmin.
        // cs[x] at scale s = min over incident edges of (r - s*0.5*D[neighbor])
        //   = inscribed slack at x once every diameter is scaled by s.
        auto scale_star = [&](std::size_t& pa, std::size_t& pb) -> double {
            double s = INF; pa = pb = 0;
            for (const FE& e : Eg) {
                const double denom = 0.5 * (D[e.i] + D[e.j]);
                if (denom <= 0.0) continue;
                const double s_e = e.r / denom;
                if (s_e < s) { s = s_e; pa = e.i; pb = e.j; }
            }
            return s;
        };
        std::vector<double> cs(N, INF);
        auto build_cs = [&](double s) {
            std::fill(cs.begin(), cs.end(), INF);
            for (const FE& e : Eg) {
                cs[e.i] = std::min(cs[e.i], e.r - s * 0.5 * D[e.j]);
                cs[e.j] = std::min(cs[e.j], e.r - s * 0.5 * D[e.i]);
            }
        };
        const std::size_t EW = event_width();
        if (swap_count_.size() != N) swap_count_.assign(N, 0);
        std::uint64_t acc = 0;
        // Cost bounds: this inner loop is O(edges) per iteration, so cap swaps per
        // pass modestly (it fires every 200 steps during the rung, so the squeeze
        // accumulates across passes) and bound rejected attempts hard.
        const long cap = std::min<long>(std::max<long>(0, p_.geo_B), 500);
        const int cap_p = std::max(1, p_.geo_max_partners);
        const long iter_cap = cap + 200;   // accepts + a small margin of rejects, then stop
        for (long iter = 0; iter < iter_cap && acc < static_cast<std::uint64_t>(cap); ++iter) {
            std::size_t a, b;
            const double s0 = scale_star(a, b);
            if (!(s0 < INF)) break;
            build_cs(s0);
            // relieve the contact pair (a,b): move the LARGER member's diameter
            // to a loose donor that accommodates it at scale s0; bring the
            // donor's smaller diameter in. Accept iff it strictly raises s*.
            bool did = false;
            for (int which = 0; which < 2 && !did; ++which) {
                const std::uint32_t c = (which == 0)
                    ? (D[a] >= D[b] ? static_cast<std::uint32_t>(a) : static_cast<std::uint32_t>(b))
                    : (D[a] >= D[b] ? static_cast<std::uint32_t>(b) : static_cast<std::uint32_t>(a));
                const double Rc = 0.5 * D[c];
                // best donor: smallest diameter among sites with room for D[c] at s0
                long best_d = -1; double best_dD = 0.0;
                const int cscan = std::max(1, cap_p);
                // scan the whole particle set with a stride cap for cost control
                const long stride = std::max<long>(1, static_cast<long>(N) / (64L * cscan));
                for (long d = 0; d < static_cast<long>(N); d += stride) {
                    if (static_cast<std::uint32_t>(d) == c) continue;
                    if (!(cs[d] < INF)) continue;            // isolated (no near edge): treat as no info
                    if (!(D[d] < D[c])) continue;            // donor must bring a SMALLER diameter
                    if (cs[d] < s0 * Rc) continue;           // site d must ACCOMMODATE D[c] at s0
                    const double dD = D[c] - D[d];
                    if (dD > best_dD) { best_dD = dD; best_d = d; }
                }
                if (best_d < 0) continue;
                const std::uint32_t d = static_cast<std::uint32_t>(best_d);
                // tentative swap, check global s* strictly increased
                std::swap(D[c], D[d]); std::swap(D0[c], D0[d]);
                std::size_t na, nb; const double s1 = scale_star(na, nb);
                if (s1 > s0 * (1.0 + 1e-12)) {
                    // keep: fix ranks + bookkeeping
                    const std::uint32_t qc = rank_of_[c], qd = rank_of_[d];
                    rank_of_[c] = qd; rank_of_[d] = qc;
                    at_rank_[qc] = d; at_rank_[qd] = c;
                    touched.push_back(c); touched.push_back(d);
                    ++swap_count_[c]; ++swap_count_[d];
                    ++acc; did = true;
                    if (p_.event_log_stride > 0 && !event_log_full_ &&
                        event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
                        const double nan = std::numeric_limits<double>::quiet_NaN();
                        event_buf_.push_back(static_cast<double>(current_step_));
                        event_buf_.push_back(static_cast<double>(c));
                        event_buf_.push_back(static_cast<double>(d));
                        event_buf_.push_back(0.0);
                        event_buf_.push_back(0.5 * (D[c] - D[d]));
                        event_buf_.push_back(s1 - s0);          // s* lift
                        event_buf_.push_back(-8.0);             // marker: H12c shrink-relief
                        event_buf_.push_back(-8.0);
                        for (std::size_t q = 8; q < EW; ++q)
                            event_buf_.push_back(q == 8 ? s0 : q == 9 ? s1 : nan);
                    }
                } else {
                    std::swap(D[c], D[d]); std::swap(D0[c], D0[d]);  // revert
                }
            }
            if (!did) break;   // no swap raised s* -> inner loop exhausted
        }
        geo_bottleneck_trades_ += acc;
        return acc;
    }
    // -----------------------------------------------------------------------

    std::vector<double> cmin(N, INF);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            cmin[ed.i] = std::min(cmin[ed.i], r - 0.5 * D[ed.j]);
            cmin[ed.j] = std::min(cmin[ed.j], r - 0.5 * D[ed.i]);
        }
    }
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    struct Cand { double gain; std::uint32_t i, j; };
    std::vector<Cand> cands;
    cands.reserve(N * 4);
    const std::uint32_t* at_rank = at_rank_.data();
    const int cap_p = std::max(1, p_.geo_max_partners);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!(cmin[i] < INF)) continue;
        const double Ri = 0.5 * D[i];
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0) continue;
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap_p - 1) / cap_p);
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!(cmin[j] < INF)) continue;
            const double Rj = 0.5 * D[j];
            if (!(Rj < Ri)) continue;                        // j must bring a SMALLER diameter to i
            if (cmin[j] < Ri) continue;                      // site j must ACCOMMODATE the bigger Ri
            // pair bottleneck before vs after swap (frozen c; swapped radii)
            const double before = std::min(cmin[i] / Ri, cmin[j] / Rj);
            const double after  = std::min(cmin[i] / Rj, cmin[j] / Ri);
            const double gain = after - before;              // > 0 raises the pair bottleneck
            if (gain <= 0.0) continue;
            cands.push_back({gain, i, j});
        }
    }
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) { return a.gain > b.gain; });
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Cand& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j]) continue;
        used[c.i] = used[c.j] = 1;
        const std::uint32_t i = c.i, j = c.j;
        const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
        const double d_rel = D0[j] / D0[i] - 1.0;
        const double dR = 0.5 * (D[j] - D[i]);
        std::swap(D[i], D[j]);
        std::swap(D0[i], D0[j]);
        rank_of_[i] = qj; rank_of_[j] = qi;
        at_rank_[qi] = j; at_rank_[qj] = i;
        touched.push_back(i); touched.push_back(j);
        ++swap_count_[i]; ++swap_count_[j];
        ++acc;
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back(d_rel);
            event_buf_.push_back(dR);
            event_buf_.push_back(c.gain);          // pair bottleneck lift
            event_buf_.push_back(-7.0);            // marker: bottleneck-relief pass
            event_buf_.push_back(-7.0);
            for (std::size_t q = 8; q < EW; ++q)
                event_buf_.push_back(q == 8 ? cmin[i] : q == 9 ? cmin[j] : nan);
        }
    }
    geo_bottleneck_trades_ += acc;
    return acc;
}

// H2 [107/117]: 3-cycle sort pass. Runs AFTER geo_pass on the same geo step.
// Objective J = sum_p D_p * cmin_p (match big D to big cage). A cyclic rotation of
// diameters among a triple (i,j,k) has dJ that is NOT the sum of its pairwise dJ, so
// it can improve a 2-exchange-stable state. We accept a rotation only when it beats
// every transposition among the triple (isolating the genuine k-opt gain). Greedy
// (T=0), each particle once, <= geo_B rotations. cmin frozen at pass start (as geo_pass).
std::uint64_t SwapEngine::geo_cycle_pass(std::size_t N, std::vector<double>& D,
                                         std::vector<double>& D0,
                                         const SwapEdgeSink& edges,
                                         std::vector<std::uint32_t>& touched)
{
    if (N < 3 || !edges.want_edges) return 0;
    if (!ranks_valid_) rebuild_ranks(D0);
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> cmin(N, INF);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            cmin[ed.i] = std::min(cmin[ed.i], r - 0.5 * D[ed.j]);
            cmin[ed.j] = std::min(cmin[ed.j], r - 0.5 * D[ed.i]);
        }
    }
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    struct Tri { double gain; std::uint32_t i, j, k; unsigned char perm; };  // perm: 0=A,1=B
    std::vector<Tri> cands;
    cands.reserve(N * 4);
    const std::uint32_t* at_rank = at_rank_.data();
    const long cap = std::max<long>(1, p_.geo_cycle_cap);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!(cmin[i] < INF)) continue;
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0 + 1) continue;                       // need at least two partners above
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap - 1) / cap);
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!(cmin[j] < INF)) continue;
            const long p0 = q + 1, p1 = static_cast<long>(ghi[q]);
            if (p1 < p0) continue;
            const long span2 = p1 - p0 + 1;
            const long stride2 = std::max<long>(1, (span2 + cap - 1) / cap);
            for (long q2 = p0; q2 <= p1; q2 += stride2) {
                const std::uint32_t k = at_rank[q2];
                if (!(cmin[k] < INF)) continue;
                // ranks ascending => D0_i <= D0_j <= D0_k; use D (=D0*kappa) for J.
                const double a = D[i], b = D[j], c = D[k];
                const double x = cmin[i], y = cmin[j], z = cmin[k];
                // rotations (positions i,j,k receive diameters):
                const double dJA = (c - a) * x + (a - b) * y + (b - c) * z;  // (c,a,b)
                const double dJB = (b - a) * x + (c - b) * y + (a - c) * z;  // (b,c,a)
                // transpositions:
                const double tij = (b - a) * (x - y);          // swap D at i,j
                const double tjk = (c - b) * (y - z);          // swap D at j,k
                const double tik = (c - a) * (x - z);          // swap D at i,k
                const double best_t = std::max(tij, std::max(tjk, tik));
                const double best_r = std::max(dJA, dJB);
                if (best_r <= 0.0 || best_r <= best_t) continue;   // cycle-only, improving
                cands.push_back({best_r, i, j, k,
                                 static_cast<unsigned char>(dJA >= dJB ? 0 : 1)});
            }
        }
    }
    std::sort(cands.begin(), cands.end(),
              [](const Tri& u, const Tri& w) { return u.gain > w.gain; });
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Tri& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j] || used[c.k]) continue;
        used[c.i] = used[c.j] = used[c.k] = 1;
        const std::uint32_t i = c.i, j = c.j, k = c.k;
        const std::uint32_t ri = rank_of_[i], rj = rank_of_[j], rk = rank_of_[k];
        const double Di = D[i], Dj = D[j], Dk = D[k];
        const double D0i = D0[i], D0j = D0[j], D0k = D0[k];
        if (c.perm == 0) {                 // perm A: positions get (Dk, Di, Dj)
            D[i] = Dk; D[j] = Di; D[k] = Dj;
            D0[i] = D0k; D0[j] = D0i; D0[k] = D0j;
            rank_of_[i] = rk; rank_of_[j] = ri; rank_of_[k] = rj;
            at_rank_[rk] = i; at_rank_[ri] = j; at_rank_[rj] = k;
        } else {                           // perm B: positions get (Dj, Dk, Di)
            D[i] = Dj; D[j] = Dk; D[k] = Di;
            D0[i] = D0j; D0[j] = D0k; D0[k] = D0i;
            rank_of_[i] = rj; rank_of_[j] = rk; rank_of_[k] = ri;
            at_rank_[rj] = i; at_rank_[rk] = j; at_rank_[ri] = k;
        }
        touched.push_back(i); touched.push_back(j); touched.push_back(k);
        ++swap_count_[i]; ++swap_count_[j]; ++swap_count_[k];
        ++acc;
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back(static_cast<double>(k));
            event_buf_.push_back(0.0);
            event_buf_.push_back(c.gain);
            event_buf_.push_back(-3.0);            // marker: 3-cycle pass
            event_buf_.push_back(-3.0);
            for (std::size_t q = 8; q < EW; ++q) event_buf_.push_back(nan);
        }
    }
    geo_cycle_trades_ += acc;
    return acc;
}

// H34 [194]: 4-cycle sort pass. Runs AFTER geo_cycle_pass on the same geo step when
// p_.geo_cycle4 != 0. Same objective J = sum_p D_p * cmin_p, cmin frozen at pass start.
// For each in-window quadruple (i<j<k<l by D0 rank) all 24 diameter permutations are
// scored; a 4-cycle is accepted only when it is improving AND strictly better than every
// lower-order permutation of the same quadruple (6 transpositions, 3 double
// transpositions, 8 3-cycles) - isolating the genuine 4-opt gain exactly as
// geo_cycle_pass isolates 3-opt from 2-opt. Greedy (T=0), each particle once,
// <= geo_B rotations. Event marker -4.0. Cost ~ N * geo_cycle4_cap^3 quadruples per pass.
namespace {
struct Perm4Table {
    unsigned char P[24][4];
    bool is4[24];
    Perm4Table() {
        unsigned char p[4] = {0, 1, 2, 3};
        int m = 0;
        do {
            bool fixed = false;
            for (int s = 0; s < 4; ++s) { P[m][s] = p[s]; if (p[s] == s) fixed = true; }
            // a 4-cycle: no fixed point and not a product of two disjoint transpositions
            is4[m] = !fixed && (p[p[0]] != 0);
            ++m;
        } while (std::next_permutation(p, p + 4) && m < 24);
    }
};
const Perm4Table& perm4_table() { static const Perm4Table t; return t; }
}  // namespace

std::uint64_t SwapEngine::geo_cycle4_pass(std::size_t N, std::vector<double>& D,
                                          std::vector<double>& D0,
                                          const SwapEdgeSink& edges,
                                          std::vector<std::uint32_t>& touched)
{
    if (N < 4 || !edges.want_edges) return 0;
    if (!ranks_valid_) rebuild_ranks(D0);
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> cmin(N, INF);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            if (r >= 1.5 * 0.5 * (D[ed.i] + D[ed.j])) continue;
            cmin[ed.i] = std::min(cmin[ed.i], r - 0.5 * D[ed.j]);
            cmin[ed.j] = std::min(cmin[ed.j], r - 0.5 * D[ed.i]);
        }
    }
    std::vector<std::uint32_t> glo, ghi;
    build_window_table_for(p_.geo_w, glo, ghi);
    const Perm4Table& PT = perm4_table();
    struct Quad { double gain; std::uint32_t i, j, k, l; unsigned char perm; };
    std::vector<Quad> cands;
    cands.reserve(N * 2);
    const std::uint32_t* at_rank = at_rank_.data();
    const long cap = std::max<long>(1, p_.geo_cycle4_cap);
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = at_rank[r];
        if (!(cmin[i] < INF)) continue;
        const long q0 = static_cast<long>(r) + 1, q1 = static_cast<long>(ghi[r]);
        if (q1 < q0 + 2) continue;                       // need at least three partners above
        const long span = q1 - q0 + 1;
        const long stride = std::max<long>(1, (span + cap - 1) / cap);
        for (long q = q0; q <= q1; q += stride) {
            const std::uint32_t j = at_rank[q];
            if (!(cmin[j] < INF)) continue;
            const long p0 = q + 1, p1 = static_cast<long>(ghi[q]);
            if (p1 < p0 + 1) continue;
            const long span2 = p1 - p0 + 1;
            const long stride2 = std::max<long>(1, (span2 + cap - 1) / cap);
            for (long q2 = p0; q2 <= p1; q2 += stride2) {
                const std::uint32_t k = at_rank[q2];
                if (!(cmin[k] < INF)) continue;
                const long s0 = q2 + 1, s1 = static_cast<long>(ghi[q2]);
                if (s1 < s0) continue;
                const long span3 = s1 - s0 + 1;
                const long stride3 = std::max<long>(1, (span3 + cap - 1) / cap);
                for (long q3 = s0; q3 <= s1; q3 += stride3) {
                    const std::uint32_t l = at_rank[q3];
                    if (!(cmin[l] < INF)) continue;
                    const double Dv[4] = {D[i], D[j], D[k], D[l]};
                    const double cv[4] = {cmin[i], cmin[j], cmin[k], cmin[l]};
                    const double base = Dv[0] * cv[0] + Dv[1] * cv[1] + Dv[2] * cv[2] + Dv[3] * cv[3];
                    double best4 = -INF, best_low = 0.0;     // identity (0) is a lower-order bound
                    int best4p = -1;
                    for (int m = 1; m < 24; ++m) {           // m=0 is the identity
                        const unsigned char* P = PT.P[m];
                        const double dJ = Dv[P[0]] * cv[0] + Dv[P[1]] * cv[1] +
                                          Dv[P[2]] * cv[2] + Dv[P[3]] * cv[3] - base;
                        if (PT.is4[m]) { if (dJ > best4) { best4 = dJ; best4p = m; } }
                        else if (dJ > best_low) best_low = dJ;
                    }
                    if (best4p < 0 || best4 <= 0.0 || best4 <= best_low) continue;   // 4-cycle-only, improving
                    cands.push_back({best4, i, j, k, l, static_cast<unsigned char>(best4p)});
                }
            }
        }
    }
    std::sort(cands.begin(), cands.end(),
              [](const Quad& u, const Quad& w) { return u.gain > w.gain; });
    std::vector<unsigned char> used(N, 0);
    std::uint64_t acc = 0;
    const std::size_t EW = event_width();
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    for (const Quad& c : cands) {
        if (acc >= static_cast<std::uint64_t>(std::max<long>(0, p_.geo_B))) break;
        if (used[c.i] || used[c.j] || used[c.k] || used[c.l]) continue;
        used[c.i] = used[c.j] = used[c.k] = used[c.l] = 1;
        const std::uint32_t idx[4] = {c.i, c.j, c.k, c.l};
        const unsigned char* P = PT.P[c.perm];
        std::uint32_t rk[4]; double dd[4], d0[4];
        for (int s = 0; s < 4; ++s) { rk[s] = rank_of_[idx[s]]; dd[s] = D[idx[s]]; d0[s] = D0[idx[s]]; }
        for (int s = 0; s < 4; ++s) {                    // slot s receives the diameter of slot P[s]
            D[idx[s]] = dd[P[s]]; D0[idx[s]] = d0[P[s]];
            rank_of_[idx[s]] = rk[P[s]]; at_rank_[rk[P[s]]] = idx[s];
            touched.push_back(idx[s]); ++swap_count_[idx[s]];
        }
        ++acc;
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + EW <= static_cast<std::size_t>(p_.event_log_cap) * EW) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(c.i));
            event_buf_.push_back(static_cast<double>(c.j));
            event_buf_.push_back(static_cast<double>(c.k));
            event_buf_.push_back(static_cast<double>(c.l));
            event_buf_.push_back(c.gain);
            event_buf_.push_back(-4.0);            // marker: 4-cycle pass
            event_buf_.push_back(-4.0);
            for (std::size_t q = 8; q < EW; ++q) event_buf_.push_back(nan);
        }
    }
    geo_cycle4_trades_ += acc;
    return acc;
}

double SwapEngine::delta_E(std::uint32_t i, std::uint32_t j) const {
    const double Ri = rad_[i], Rj = rad_[j];
    if (Ri == Rj) return 0.0;
    double dE = 0.0;
    // particle a takes radius Rnew (was Rold); the a-partner contact is unchanged.
    auto side = [&](std::uint32_t a, double Rold, double Rnew, std::uint32_t partner) {
        const std::uint32_t lo = adj_off_[a], hi = adj_off_[a + 1];
        double s = 0.0;
        for (std::uint32_t q = lo; q < hi; ++q) {
            const std::uint32_t k = adj_idx_[q];
            if (k == partner) continue;
            const double r = adj_r_[q];
            const double Rk = rad_[k];
            s += e_term(Rnew + Rk - r) - e_term(Rold + Rk - r);
        }
        return s;
    };
    dE += side(i, Ri, Rj, j);
    dE += side(j, Rj, Ri, i);
    return dE;
}

double SwapEngine::delta_E3(std::uint32_t i, std::uint32_t j, std::uint32_t k) const {
    // rotation: i takes R_j, j takes R_k, k takes R_i. Unlike a pair swap, the
    // intra-trio pair sums change, so trio-internal contacts are NOT skipped: each
    // is counted once with new radii on both sides.
    const std::uint32_t trio[3] = {i, j, k};
    const double Rold[3] = {rad_[i], rad_[j], rad_[k]};
    const double Rnew[3] = {rad_[j], rad_[k], rad_[i]};
    auto newR = [&](std::uint32_t a, bool& is_trio) -> double {
        for (int t = 0; t < 3; ++t) if (trio[t] == a) { is_trio = true; return Rnew[t]; }
        is_trio = false; return rad_[a];
    };
    double dE = 0.0;
    for (int t = 0; t < 3; ++t) {
        const std::uint32_t a = trio[t];
        const std::uint32_t lo = adj_off_[a], hi = adj_off_[a + 1];
        for (std::uint32_t q = lo; q < hi; ++q) {
            const std::uint32_t b = adj_idx_[q];
            const double r = adj_r_[q];
            bool b_trio = false;
            const double Rb_new = newR(b, b_trio);
            if (b_trio) {
                if (b <= a) continue;   // count each intra-trio edge once
                // both sides change
                double Rb_old = rad_[b];
                dE += e_term(Rnew[t] + Rb_new - r) - e_term(Rold[t] + Rb_old - r);
            } else {
                dE += e_term(Rnew[t] + rad_[b] - r) - e_term(Rold[t] + rad_[b] - r);
            }
        }
    }
    return dE;
}

double SwapEngine::edge_margin(const std::vector<double>& D, const std::vector<double>& D0) {
    const std::size_t N = D.size();
    if (!ranks_valid_) rebuild_ranks(D0);
    double delta_max = 0.0;
    if (p_.mode == 1) {
        const std::size_t W = static_cast<std::size_t>(std::max(1, p_.rank_window));
        for (std::size_t r = 0; r + W < N; ++r) {
            const double dd = 0.5 * (D[at_rank_[r + W]] - D[at_rank_[r]]);
            if (dd > delta_max) delta_max = dd;
        }
        if (W >= N && N > 0) delta_max = 0.5 * (D[at_rank_[N - 1]] - D[at_rank_[0]]);
    } else if (p_.mode == 5) {
        double dmax = 0.0;
        for (std::size_t i = 0; i < N; ++i) if (D[i] > dmax) dmax = D[i];
        // ctrl mode: the adaptive window can sit above window_rel (up to w_max);
        // cover whichever is in effect so the collected edges stay sufficient
        double w_eff = p_.window_rel;
        if (p_.ctrl && w_cur_ > w_eff) w_eff = w_cur_;
        delta_max = 0.5 * w_eff * dmax;
        // rescue channel: one extra move of up to rescue_w per particle per pass;
        // margin is only requested on rescue steps, so this widening is paid only
        // there (the cap factor below covers the main channel; add the single
        // rescue move on top after the cap multiplication via delta_rescue_)
        if ((p_.rescue > 0.0 || p_.cycle3 > 0.0) && p_.rescue_w > 0.0) {
            const int cap_ = std::max(1, p_.max_per_particle);
            // fold the one rescue move into delta_max pre-cap: cap*delta_main +
            // 1*delta_rescue == cap*(delta_main + delta_rescue/cap)
            delta_max += 0.5 * p_.rescue_w * dmax / static_cast<double>(cap_);
        }
    } else {
        double meanD = 0.0;
        for (std::size_t i = 0; i < N; ++i) meanD += D[i];
        meanD /= static_cast<double>(std::max<std::size_t>(1, N));
        delta_max = 0.5 * p_.size_tol * meanD;
    }
    // A particle may be swapped at most `cap` times per sweep, so its radius moves by
    // at most cap*delta_max within the sweep; using that as the collection margin
    // keeps every sequential dE exact.
    const int cap = std::max(1, p_.max_per_particle);
    return static_cast<double>(cap) * delta_max;
}

// ---------------------------------------------------------------------------
// Proxy sweep: O(1) per candidate. With g_a = sum of a's absolute overlaps and z_a
// its contact count (both from the force loop), a swap where `grow` gains dR and
// `shrink` loses dR changes the objective by approximately
//     dE ~ (g_grow - g_shrink) dR + 1/2 (z_grow + z_shrink) dR^2
// (exact if no contact is created or lost during the change). Biased towards
// acceptance: created contacts and lost contacts both make the true dE larger.
// After an accepted swap the pair's own (g, z) are propagated to first order.
SwapSweepStats SwapEngine::sweep_proxy(std::size_t N,
                                       std::vector<double>& D,
                                       std::vector<double>& D0,
                                       const double* gsum_in,
                                       const std::vector<std::size_t>& z,
                                       std::vector<std::uint32_t>& touched)
{
    SwapSweepStats st;
    const auto t0 = std::chrono::steady_clock::now();
    rad_.resize(N);
    gsum_.resize(N);
    zc_.resize(N);
    double gmean = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        rad_[i] = 0.5 * D[i];
        gsum_[i] = gsum_in[i];
        zc_[i] = static_cast<std::uint32_t>(z[i]);
        gmean += gsum_in[i];
    }
    gmean /= static_cast<double>(N);
    st.contacts = 0;
    for (std::size_t i = 0; i < N; ++i) st.contacts += z[i];
    st.contacts /= 2;
    st.e_scale = 0.0;
    per_particle_count_.assign(N, 0u);
    const std::uint32_t cap = static_cast<std::uint32_t>(std::max(1, p_.max_per_particle));
    const auto t1 = std::chrono::steady_clock::now();

    const std::uint64_t n_cand = static_cast<std::uint64_t>(p_.swap * static_cast<double>(N) + 0.5);
    const std::uint32_t* rank_of = rank_of_.data();
    const std::uint32_t* at_rank = at_rank_.data();
    double* rad = rad_.data();
    double* gsum = gsum_.data();
    std::uint32_t* zc = zc_.data();
    std::uint32_t* ppc = per_particle_count_.data();
    const std::uint64_t W = static_cast<std::uint64_t>(std::max(1, p_.rank_window));
    const std::uint64_t phase = (p_.mode == 3) ? (next_u64() % N) : 0;
    for (std::uint64_t c = 0; c < n_cand; ++c) {
        const std::uint64_t u = next_u64();
        std::uint32_t i, j;
        long ri;
        if (p_.mode == 3) { ri = static_cast<long>((c + phase) % N); i = at_rank[ri]; }
        else { i = static_cast<std::uint32_t>(u % N); ri = static_cast<long>(rank_of[i]); }
        long r;
        if (p_.mode == 5) {
            // directed, scale-free window: partner rank uniform on the side chosen by
            // pressure within [win_lo, win_hi] of rank ri
            const long lo = win_lo_[ri], hi = win_hi_[ri];
            const bool shrink = gsum[i] > gmean;
            const long a = shrink ? lo : ri + 1, b = shrink ? ri - 1 : hi;
            ++st.candidates;
            if (b < a) continue;
            r = a + static_cast<long>((u >> 20) % static_cast<std::uint64_t>(b - a + 1));
        } else {
            long off = 1 + static_cast<long>((u >> 20) % W);
            if (p_.mode == 4) { if (gsum[i] > gmean) off = -off; }
            else if (u >> 63) off = -off;
            r = ri + off;
            ++st.candidates;
            if (r < 0 || r >= static_cast<long>(N)) continue;
        }
        j = at_rank[static_cast<std::size_t>(r)];
        if (i == j) continue;
        if (ppc[i] >= cap || ppc[j] >= cap) continue;
        const bool i_grows = rad[i] < rad[j];
        const std::uint32_t a = i_grows ? i : j;   // grows
        const std::uint32_t b = i_grows ? j : i;   // shrinks
        const double dR = rad[b] - rad[a];
        if (dR <= 0.0) continue;
        ++st.evaluated;
        const double dE = (gsum[a] - gsum[b]) * dR + 0.5 * static_cast<double>(zc[a] + zc[b]) * dR * dR;
        if (dE < 0.0) {
            std::swap(rad[i], rad[j]);
            std::swap(D[i], D[j]);
            std::swap(D0[i], D0[j]);
            const std::uint32_t rri = rank_of_[i], rrj = rank_of_[j];
            rank_of_[i] = rrj; rank_of_[j] = rri;
            at_rank_[rri] = j; at_rank_[rrj] = i;
            touched.push_back(i);
            touched.push_back(j);
            ++ppc[i]; ++ppc[j];
            if (swap_count_.size() != N) swap_count_.assign(N, 0);
            ++swap_count_[i]; ++swap_count_[j];
            ++st.accepted;
            st.dE_sum += dE;
            // first-order propagation of the pair's own overlap sums
            gsum[a] += static_cast<double>(zc[a]) * dR;
            gsum[b] = std::max(0.0, gsum[b] - static_cast<double>(zc[b]) * dR);
        }
    }
    const auto t2 = std::chrono::steady_clock::now();
    t_build_ += std::chrono::duration<double>(t1 - t0).count();
    t_cand_ += std::chrono::duration<double>(t2 - t1).count();
    n_eval_total_ += st.evaluated; n_cand_total_ += st.candidates; ++n_sweeps_;
    return st;
}

// ---------------------------------------------------------------------------
// Rank-space proxy sweep. Same acceptance rule as sweep_proxy, but organised so
// that all per-candidate work is streaming in rank space:
//   1. scatter g, z of every particle to its rank (one pass, N random stores);
//   2. walk the ranks in a strided order; the partner rank is a few positions away,
//      so its data sits in the same cache lines; on acceptance swap the two
//      particles' rank entries and their D0/D (the only random writes).
// Per candidate this costs O(1) with ~no dependent cache misses at any N, which is
// what keeps the sweep's cost per particle constant when N no longer fits in cache.
SwapSweepStats SwapEngine::sweep_proxy_rank(std::size_t N,
                                            std::vector<double>& D,
                                            std::vector<double>& D0,
                                            const double* gsum_in,
                                            const std::vector<std::size_t>& z,
                                            std::vector<std::uint32_t>& touched)
{
    // Rank-banded, threaded sweep (speed-first; no cross-thread bit contract).
    // Rank space is partitioned circularly into bands whose boundaries shift by a
    // random offset each sweep; every candidate's partner window is clipped to its
    // band, so a band touches only its own rank entries and the particles that own
    // them -- bands are independent jobs (no locks, no shared writes) scheduled
    // dynamically over the OpenMP threads. Each band draws from its own splitmix
    // stream seeded (run stream, sweep, band). Boundary clipping loses a fraction
    // ~W/B of window area per sweep; the rotation redistributes it across ranks.
    SwapSweepStats st;
    const auto t0 = std::chrono::steady_clock::now();
    const double kappa = (D0[0] > 0.0) ? D[0] / D0[0] : 1.0;   // D = D0 * kappa (uniform)
    rank_g_.resize(N); rank_z_.resize(N); rank_ppc_.resize(N);
    double gmean = 0.0;
    std::uint64_t zsum = 0;
    {
        // rank_of_ is a permutation: each rank written exactly once -> no conflicts
        double* rg0 = rank_g_.data(); std::uint32_t* rz0 = rank_z_.data();
        std::uint32_t* rp0 = rank_ppc_.data();
        const std::uint32_t* ro = rank_of_.data();
        #pragma omp parallel for reduction(+:gmean, zsum) schedule(static)
        for (std::size_t i = 0; i < N; ++i) {
            const std::uint32_t r = ro[i];
            rg0[r] = gsum_in[i];
            rz0[r] = static_cast<std::uint32_t>(z[i]);
            rp0[r] = 0u;
            gmean += gsum_in[i];
            zsum += z[i];
        }
    }
    gmean /= static_cast<double>(N);
    st.contacts = static_cast<std::size_t>(zsum / 2);
    if (swap_count_.size() != N) swap_count_.assign(N, 0);
    const std::uint32_t cap = static_cast<std::uint32_t>(std::max(1, p_.max_per_particle));

    // ---- band layout ---------------------------------------------------------
    // Band width: at least twice the largest window span (so clipping is a
    // boundary effect, not the rule) and large enough that bands stay meaty;
    // band count targets ~4 jobs per thread for dynamic balance.
    std::size_t wmax = wmax_cached_;
    if (wmax == 0) {
        wmax = 1;
        if (p_.mode == 5 && !win_lo_.empty()) {
            for (std::size_t r = 0; r < N; ++r) {
                const std::size_t span = static_cast<std::size_t>(win_hi_[r] - win_lo_[r] + 1);
                if (span > wmax) wmax = span;
            }
        } else {
            wmax = static_cast<std::size_t>(std::max(1, p_.rank_window)) * 2 + 1;
        }
        wmax_cached_ = wmax;
    }
    const int nthreads = omp_get_max_threads();
    std::size_t nb_target = static_cast<std::size_t>(std::max(1, 4 * nthreads));
    std::size_t B = std::max<std::size_t>(2 * wmax, (N + nb_target - 1) / nb_target);
    if (B > N) B = N;
    const std::size_t nbands = (N + B - 1) / B;
    const std::uint64_t sweep_salt = next_u64();
    const std::uint64_t shift = (nbands > 1) ? (sweep_salt % B) : 0;

    const std::uint32_t* at0 = at_rank_.data();
    std::uint32_t* at_rank = at_rank_.data();
    std::uint32_t* rank_of = rank_of_.data();
    const double* rD0 = rank_D0_.data();
    double* rg = rank_g_.data();
    std::uint32_t* rz = rank_z_.data();
    std::uint32_t* rppc = rank_ppc_.data();
    const std::uint64_t W = static_cast<std::uint64_t>(std::max(1, p_.rank_window));
    (void)at0;

    std::uint64_t tot_cand = 0, tot_eval = 0, tot_acc = 0;
    double tot_dE = 0.0;
    std::vector<std::vector<std::uint32_t>> touched_tls(nbands);
    std::vector<std::vector<double>> events_tls(
        p_.event_log_stride > 0 ? nbands : 0);
    if (p_.event_log_stride <= 0) events_tls.resize(nbands);  // keep indexing valid

    #pragma omp parallel for schedule(dynamic) \
        reduction(+:tot_cand, tot_eval, tot_acc, tot_dE)
    for (std::size_t b = 0; b < nbands; ++b) {
        // circular band [lo, hi): ranks lo..hi-1, where the partition is offset by
        // `shift`. Band b covers [shift + b*B, shift + (b+1)*B) mod N -- to keep the
        // arithmetic simple we clamp to the linear range and give the wrap segment
        // [0, shift) to the last band.
        std::size_t lo = shift + b * B;
        std::size_t hi = std::min(N, lo + B);
        const bool last = (b == nbands - 1);
        // wrap segment handled by the last band as a second range [0, shift)
        const std::size_t lo2 = 0, hi2 = last ? shift : 0;
        if (lo >= N) { lo = 0; hi = 0; }
        const std::size_t size1 = hi - lo, size2 = hi2 - lo2;
        const std::size_t bsize = size1 + size2;
        if (bsize < 2) continue;
        std::vector<std::uint32_t>& tch = touched_tls[b];
        std::vector<double>& ev = events_tls[b];
        std::uint64_t ev_ctr = 0;
        // per-band RNG: splitmix64 on a per-band state
        std::uint64_t rs = sweep_salt ^ (0x9E3779B97F4A7C15ULL * (b + 1));
        auto rnd = [&rs]() {
            std::uint64_t x = (rs += 0x9E3779B97F4A7C15ULL);
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
            return x ^ (x >> 31);
        };
        // strided walk over the band's ranks
        std::uint64_t stride = (static_cast<std::uint64_t>(std::sqrt(double(bsize))) | 1u);
        auto gcd = [](std::uint64_t a, std::uint64_t c) { while (c) { const std::uint64_t t = a % c; a = c; c = t; } return a; };
        while (gcd(stride, bsize) != 1) stride += 2;
        std::uint64_t walk = rnd() % bsize;
        auto rank_at = [&](std::uint64_t k) -> std::size_t {
            return (k < size1) ? (lo + k) : (lo2 + (k - size1));
        };
        auto in_band = [&](long r) -> bool {
            if (r < 0 || r >= static_cast<long>(N)) return false;
            const std::size_t ur = static_cast<std::size_t>(r);
            return (ur >= lo && ur < hi) || (ur >= lo2 && ur < hi2);
        };
        // draw budget: fixed s (default) or ctrl quota (draw until a*bsize accepted
        // or s_max*bsize drawn). With the window pinned at its floor the accept
        // quota is unreachable (acceptance is physics-limited, far below target_R)
        // and drawing to the cap buys ~4x the draws of a healthy sweep for a small
        // fraction of the exchanges [log 031]; floor_budget caps such sweeps at the
        // healthy-phase budget (target_a/target_R per particle) instead.
        const bool at_floor = p_.ctrl && w_cur_ >= 0.0 &&
                              (w_cur_ <= p_.w_min * (1.0 + 1e-9));
        const double ctrl_pp = (at_floor && p_.floor_budget > 0.0)
            ? std::min(p_.s_max, p_.floor_budget * p_.target_a
                                     / std::max(p_.target_R, 1e-9))
            : p_.s_max;
        const std::uint64_t n_cand = p_.ctrl
            ? static_cast<std::uint64_t>(ctrl_pp * double(bsize) + 0.5)
            : static_cast<std::uint64_t>(p_.swap * double(bsize) + 0.5);
        const std::uint64_t quota = p_.ctrl
            ? static_cast<std::uint64_t>(p_.target_a * double(bsize) + 0.5)
            : ~0ull;
        std::uint64_t cand = 0, eval = 0, acc = 0; double dEs = 0.0;
        for (std::uint64_t c = 0; c < n_cand && acc < quota; ++c) {
            const std::uint64_t u = rnd();
            const long ri = static_cast<long>(rank_at(walk));
            walk += stride; if (walk >= bsize) walk -= bsize;
            ++cand;
            const bool shrink = (p_.mode == 5 || p_.mode == 4) ? (rg[ri] > gmean) : ((u >> 63) != 0);
            long rj;
            if (p_.mode == 5) {
                long a = shrink ? static_cast<long>(win_lo_[ri]) : ri + 1;
                long bnd = shrink ? ri - 1 : static_cast<long>(win_hi_[ri]);
                // clip to band (both segments are contiguous; the window around ri
                // stays within ri's own segment except at segment edges)
                if (static_cast<std::size_t>(ri) >= lo && static_cast<std::size_t>(ri) < hi) {
                    a = std::max(a, static_cast<long>(lo)); bnd = std::min(bnd, static_cast<long>(hi) - 1);
                } else {
                    a = std::max(a, static_cast<long>(lo2)); bnd = std::min(bnd, static_cast<long>(hi2) - 1);
                }
                if (bnd < a) continue;
                rj = a + static_cast<long>((u >> 20) % static_cast<std::uint64_t>(bnd - a + 1));
            } else {
                long off = 1 + static_cast<long>((u >> 20) % W);
                rj = shrink ? ri - off : ri + off;
                if (!in_band(rj)) continue;
            }
            if (rj == ri) continue;
            if (rppc[ri] >= cap || rppc[rj] >= cap) continue;
            const long ra = ri < rj ? ri : rj;   // grows
            const long rb = ri < rj ? rj : ri;   // shrinks
            const double dR = 0.5 * kappa * (rD0[rb] - rD0[ra]);
            if (dR <= 0.0) continue;
            ++eval;
            const double dE = (rg[ra] - rg[rb]) * dR + 0.5 * static_cast<double>(rz[ra] + rz[rb]) * dR * dR;
            if (dE < 0.0) {
                const std::uint32_t i = at_rank[ri], j = at_rank[rj];
                if (p_.event_log_stride > 0 && !event_log_full_ &&
                    (ev_ctr++ % static_cast<std::uint64_t>(p_.event_log_stride)) == 0) {
                    ev.push_back(static_cast<double>(current_step_));
                    ev.push_back(static_cast<double>(at_rank[ra]));   // grows  [067: role-ordered]
                    ev.push_back(static_cast<double>(at_rank[rb]));   // shrinks
                    ev.push_back(rD0[rb] / rD0[ra] - 1.0);   // delta_rel
                    ev.push_back(dR);                        // delta_abs (radius)
                    ev.push_back(-dE);                       // relief (proxy)
                    ev.push_back(static_cast<double>(rz[ra]));
                    ev.push_back(static_cast<double>(rz[rb]));
                    if (p_.event_log_geom) {                  // ME0-iv
                        const std::uint32_t ga = at_rank[ra], gb = at_rank[rb];
                        const bool have = cage_.size() == N;
                        const double nan = std::numeric_limits<double>::quiet_NaN();
                        ev.push_back(have ? cage_[ga] : nan);
                        ev.push_back(have ? cage_[gb] : nan);
                        ev.push_back(have ? fact_[ga] : nan);
                        ev.push_back(have ? fact_[gb] : nan);
                    }
                }
                std::swap(D[i], D[j]);
                std::swap(D0[i], D0[j]);
                at_rank[ri] = j; at_rank[rj] = i;
                rank_of[i] = static_cast<std::uint32_t>(rj); rank_of[j] = static_cast<std::uint32_t>(ri);
                const double ga = rg[ra] + static_cast<double>(rz[ra]) * dR;
                const double gb = std::max(0.0, rg[rb] - static_cast<double>(rz[rb]) * dR);
                rg[ra] = gb; rg[rb] = ga;
                std::swap(rz[ra], rz[rb]);
                const std::uint32_t pa = rppc[ra] + 1, pb = rppc[rb] + 1;
                rppc[ra] = pb; rppc[rb] = pa;
                tch.push_back(i); tch.push_back(j);
                ++swap_count_[i]; ++swap_count_[j];
                ++acc; dEs += dE;
            }
        }
        tot_cand += cand; tot_eval += eval; tot_acc += acc; tot_dE += dEs;
    }
    st.candidates = tot_cand; st.evaluated = tot_eval; st.accepted = tot_acc; st.dE_sum = tot_dE;
    for (auto& v : touched_tls) touched.insert(touched.end(), v.begin(), v.end());
    if (p_.event_log_stride > 0 && !event_log_full_) {
        for (auto& v : events_tls) {
            if (event_buf_.size() + v.size() >
                static_cast<std::size_t>(p_.event_log_cap) * event_width()) {
                event_log_full_ = true;
                break;
            }
            event_buf_.insert(event_buf_.end(), v.begin(), v.end());
        }
    }
    // ---- (a, R) controller update -------------------------------------------
    if (p_.ctrl) {
        st.w_used = w_cur_;
        const double A = tot_cand ? double(tot_acc) / double(tot_cand) : 0.0;
        // Starvation = the channel is genuinely dry, not merely under quota: the
        // window is at its floor AND the sweep accepted less than the same idle
        // fraction the fixed-window mode uses. (Quota-unmet-at-floor alone parks the
        // sweeps mid-compression at large S, where the reachable acceptance plateau
        // can sit below any fixed R target -- measured at S=20, SURF study.)
        st.starved = (tot_acc < static_cast<std::uint64_t>(p_.idle_frac * double(N))) &&
                     (w_cur_ <= p_.w_min * (1.0 + 1e-9));
        double w_new = w_cur_ * std::exp(p_.ctrl_eta * (A - p_.target_R));
        if (w_new < p_.w_min) w_new = p_.w_min;
        if (w_new > p_.w_max) w_new = p_.w_max;
        // rebuild the window table only on a material change (>2%) — O(N) two-pointer
        if (std::fabs(w_new / w_cur_ - 1.0) > 0.02) {
            w_cur_ = w_new;
            rebuild_window_table(w_cur_);
        } else {
            w_cur_ = w_new;
        }
    }
    const auto t2 = std::chrono::steady_clock::now();
    t_build_ += 0.0;
    t_cand_ += std::chrono::duration<double>(t2 - t0).count();
    n_eval_total_ += st.evaluated; n_cand_total_ += st.candidates; ++n_sweeps_;
    return st;
}

SwapSweepStats SwapEngine::sweep(std::size_t N,
                                 std::vector<double>& D,
                                 std::vector<double>& D0,
                                 const SwapEdgeSink& edges,
                                 const std::vector<std::size_t>& z,
                                 std::vector<std::uint32_t>& touched)
{
    SwapSweepStats st;
    if (N < 2 || p_.swap <= 0.0) return st;
    if (!ranks_valid_) rebuild_ranks(D0);
    if (p_.event_log_geom && p_.event_log_stride > 0) {
        if (edges.want_edges) geom_from_edges(N, D, edges, z);
        else { cage_.clear(); fact_.clear(); }
    }
    if (!p_.exact) {
        SwapSweepStats st_main;
        if (p_.ctrl || (p_.rank_space && (p_.mode == 5 || p_.mode == 4 || p_.mode == 1 || p_.mode == 3)))
            st_main = sweep_proxy_rank(N, D, D0, edges.gsum.data(), z, touched);
        else
            st_main = sweep_proxy(N, D, D0, edges.gsum.data(), z, touched);
        // rescue channel: exact-dE large-window pass; runs only on steps where the
        // force loop collected edges (rcp_generator gates that on rescue_interval)
        if ((p_.rescue > 0.0 || p_.cycle3 > 0.0) && p_.mode == 5 && edges.want_edges)
            rescue_pass(N, D, D0, edges, touched, st_main);
        return st_main;
    }
    const auto t0 = std::chrono::steady_clock::now();

    rad_.resize(N);
    for (std::size_t i = 0; i < N; ++i) rad_[i] = 0.5 * D[i];
    double meanD = 0.0;
    for (std::size_t i = 0; i < N; ++i) meanD += D[i];
    meanD /= static_cast<double>(N);
    const double dtol = p_.size_tol * meanD;
    const std::uint32_t cap = static_cast<std::uint32_t>(std::max(1, p_.max_per_particle));

    // ---- CSR adjacency over the near-pair edges (both directions) --------------
    deg_.assign(N + 1, 0u);
    std::size_t n_edges = 0;
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        n_edges += n;
        for (std::size_t q = 0; q < n; ++q) { ++deg_[v[q].i]; ++deg_[v[q].j]; }
    }
    adj_off_.assign(N + 1, 0u);
    for (std::size_t i = 0; i < N; ++i) adj_off_[i + 1] = adj_off_[i] + deg_[i];
    const std::uint32_t total = adj_off_[N];
    adj_idx_.resize(total);
    adj_r_.resize(total);
    fill_.assign(N, 0u);
    gsum_.assign(N, 0.0);
    zc_.assign(N, 0u);
    double e = 0.0; std::size_t nc = 0;
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            const std::uint32_t pi = adj_off_[ed.i] + fill_[ed.i]++;
            const std::uint32_t pj = adj_off_[ed.j] + fill_[ed.j]++;
            adj_idx_[pi] = ed.j; adj_r_[pi] = r;
            adj_idx_[pj] = ed.i; adj_r_[pj] = r;
            const double o = rad_[ed.i] + rad_[ed.j] - r;
            if (o > 0.0) { e += e_term(o); ++nc; gsum_[ed.i] += o; gsum_[ed.j] += o; ++zc_[ed.i]; ++zc_[ed.j]; }
        }
    }
    // mean overlap sum per particle: the directed mode's grow/shrink pivot
    double gmean = 0.0;
    for (std::size_t i = 0; i < N; ++i) gmean += gsum_[i];
    gmean /= static_cast<double>(N);
    (void)n_edges;
    st.contacts = nc;
    st.e_scale = nc ? e / static_cast<double>(nc) : 0.0;
    const double thresh = -p_.eps * st.e_scale;

    per_particle_count_.assign(N, 0u);
    const auto t1 = std::chrono::steady_clock::now();

    // ---- candidate loop --------------------------------------------------------
    // mode 1: i uniform, partner = rank(i) +- off, off uniform in [1, W].
    // mode 3: systematic rank sweep -- k-th candidate uses rank (k + phase) mod N as i
    //         (every rank once per N candidates), partner as in mode 1. Fewer wasted
    //         draws and sequential access to the rank table.
    const std::uint64_t n_cand = static_cast<std::uint64_t>(p_.swap * static_cast<double>(N) + 0.5);
    const std::uint32_t* rank_of = rank_of_.data();
    const std::uint32_t* at_rank = at_rank_.data();
    double* rad = rad_.data();
    double* gsum = gsum_.data();
    const std::uint32_t* zc = zc_.data();
    std::uint32_t* ppc = per_particle_count_.data();
    const std::uint64_t W = static_cast<std::uint64_t>(std::max(1, p_.rank_window));
    const std::uint64_t phase = (p_.mode == 3) ? (next_u64() % N) : 0;
    for (std::uint64_t c = 0; c < n_cand; ++c) {
        const std::uint64_t u = next_u64();          // one draw per candidate
        std::uint32_t i, j;
        if (p_.mode == 1 || p_.mode == 3 || p_.mode == 4 || p_.mode == 5) {
            long ri;
            if (p_.mode == 3) {
                ri = static_cast<long>((c + phase) % N);
                i = at_rank[ri];
            } else {
                i = static_cast<std::uint32_t>(u % N);
                ri = static_cast<long>(rank_of[i]);
            }
            long r;
            if (p_.mode == 5) {
                const long lo = win_lo_[ri], hi = win_hi_[ri];
                const double gi = gsum[i];
                const bool shrink = (gi < 0.0) ? ((u >> 63) != 0) : (gi > gmean);
                const long a = shrink ? lo : ri + 1, b = shrink ? ri - 1 : hi;
                ++st.candidates;
                if (b < a) continue;
                r = a + static_cast<long>((u >> 20) % static_cast<std::uint64_t>(b - a + 1));
            } else {
                long off = 1 + static_cast<long>((u >> 20) % W);
                if (p_.mode == 4) {
                    // directed: over-pressured particles propose to shrink (lower rank),
                    // under-pressured ones to grow; unknown pressure -> random direction
                    const double gi = gsum[i];
                    if (gi < 0.0) { if (u >> 63) off = -off; }
                    else if (gi > gmean) off = -off;
                } else {
                    if (u >> 63) off = -off;
                }
                r = ri + off;
                ++st.candidates;
                if (r < 0 || r >= static_cast<long>(N)) continue;
            }
            j = at_rank[static_cast<std::size_t>(r)];
        } else if (p_.mode == 2) {
            i = static_cast<std::uint32_t>(u % N);
            const std::uint32_t lo = adj_off_[i], hi = adj_off_[i + 1];
            ++st.candidates;
            if (hi == lo) continue;
            j = adj_idx_[lo + static_cast<std::uint32_t>((u >> 20) % (hi - lo))];
        } else {
            i = static_cast<std::uint32_t>(u % N);
            j = static_cast<std::uint32_t>((u >> 20) % N);
            ++st.candidates;
        }
        if (i == j) continue;
        if (p_.mode == 0 || p_.mode == 2) { if (std::fabs(D[i] - D[j]) > dtol) continue; }
        if (ppc[i] >= cap || ppc[j] >= cap) continue;
        // O(1) necessary condition. With g_a = sum of a's current overlaps, the exact
        // change obeys dE >= (g_grow - g_shrink) * dR  (dR > 0: convexity of the
        // growing side, and losing contacts only makes the shrinking side's decrease
        // smaller). So a swap can lower E only if the particle that SHRINKS carries
        // more overlap than the one that grows. Never rejects a dE < 0 move; the exact
        // dE below still decides. gsum < 0 marks "unknown" (near an earlier swap this
        // sweep) -> no filtering.
        {
            const bool i_grows = rad[i] < rad[j];
            const double g_grow = i_grows ? gsum[i] : gsum[j];
            const double g_shrink = i_grows ? gsum[j] : gsum[i];
            // + the growing side's second-order term over its existing contacts:
            // dE >= (g_grow - g_shrink) dR + 1/2 z_grow dR^2   (new contacts only add)
            const double dR = std::fabs(rad[i] - rad[j]);
            const double z_grow = static_cast<double>(i_grows ? zc[i] : zc[j]);
            if (p_.eps >= 0.0 && g_shrink >= 0.0 && g_grow >= 0.0 &&
                g_shrink <= g_grow + 0.5 * z_grow * dR) continue;
        }
        ++st.evaluated;
        const double dE = delta_E(i, j);
        if (dE < thresh) {   // eps > 0 stricter descent; eps < 0 admits bounded uphill moves
            std::swap(rad[i], rad[j]);
            std::swap(D[i], D[j]);
            std::swap(D0[i], D0[j]);
            const std::uint32_t ri = rank_of_[i], rj = rank_of_[j];
            rank_of_[i] = rj; rank_of_[j] = ri;
            at_rank_[ri] = j; at_rank_[rj] = i;
            touched.push_back(i);
            touched.push_back(j);
            ++ppc[i]; ++ppc[j];
            if (swap_count_.size() != N) swap_count_.assign(N, 0);
            ++swap_count_[i]; ++swap_count_[j];
            ++st.accepted;
            st.dE_sum += dE;
            // gsum of i, j and their near neighbours is now stale. Refresh i, j
            // exactly (cheap, near-lists) and mark their neighbours "unknown" (-1),
            // which the prefilter treats as "cannot exclude" -> exact dE decides.
            for (std::uint32_t q = adj_off_[i]; q < adj_off_[i + 1]; ++q) gsum[adj_idx_[q]] = -1.0;
            for (std::uint32_t q = adj_off_[j]; q < adj_off_[j + 1]; ++q) gsum[adj_idx_[q]] = -1.0;
            for (std::uint32_t a : {i, j}) {
                double g = 0.0; std::uint32_t z = 0;
                for (std::uint32_t q = adj_off_[a]; q < adj_off_[a + 1]; ++q) {
                    const double o = rad[a] + rad[adj_idx_[q]] - adj_r_[q];
                    if (o > 0.0) { g += o; ++z; }
                }
                gsum[a] = g; zc_[a] = z;
            }
        }
    }
    const auto t2 = std::chrono::steady_clock::now();
    t_build_ += std::chrono::duration<double>(t1 - t0).count();
    t_cand_ += std::chrono::duration<double>(t2 - t1).count();
    n_edges_total_ += n_edges; n_eval_total_ += st.evaluated; n_cand_total_ += st.candidates; ++n_sweeps_;
    return st;
}

// ---------------------------------------------------------------------------
// Rescue channel [log 034]: a rare SERIAL pass after the main proxy sweep, on
// steps where the force loop collected near-pair edges. Candidates come from a
// wide relative-D window (rescue_w >> the main window) and are decided by the
// EXACT dE over the collected edges; the O(1) proxy's verdict on the same
// candidate is recorded, so the pass doubles as a direct measurement of the
// proxy's false-reject rate on large-delta moves. Runs at fixed positions,
// after the parallel region; mutates D/D0 like any accepted exchange.
void SwapEngine::rescue_pass(std::size_t N, std::vector<double>& D,
                             std::vector<double>& D0, const SwapEdgeSink& edges,
                             std::vector<std::uint32_t>& touched,
                             SwapSweepStats& st)
{
    // ---- CSR adjacency + current per-particle overlap state over the edges ----
    rad_.resize(N);
    for (std::size_t i = 0; i < N; ++i) rad_[i] = 0.5 * D[i];
    deg_.assign(N + 1, 0u);
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) { ++deg_[v[q].i]; ++deg_[v[q].j]; }
    }
    adj_off_.assign(N + 1, 0u);
    for (std::size_t i = 0; i < N; ++i) adj_off_[i + 1] = adj_off_[i] + deg_[i];
    adj_idx_.resize(adj_off_[N]);
    adj_r_.resize(adj_off_[N]);
    fill_.assign(N, 0u);
    gsum_.assign(N, 0.0);
    zc_.assign(N, 0u);
    double e = 0.0; std::size_t nc = 0;
    for (std::size_t t = 0; t < edges.per_thread.size(); ++t) {
        const SwapEdge* v = edges.per_thread[t].data();
        const std::size_t n = edges.count[t];
        for (std::size_t q = 0; q < n; ++q) {
            const SwapEdge& ed = v[q];
            const double r = std::sqrt(ed.d2);
            const std::uint32_t pi = adj_off_[ed.i] + fill_[ed.i]++;
            const std::uint32_t pj = adj_off_[ed.j] + fill_[ed.j]++;
            adj_idx_[pi] = ed.j; adj_r_[pi] = r;
            adj_idx_[pj] = ed.i; adj_r_[pj] = r;
            const double o = rad_[ed.i] + rad_[ed.j] - r;
            if (o > 0.0) { e += e_term(o); ++nc; gsum_[ed.i] += o; gsum_[ed.j] += o; ++zc_[ed.i]; ++zc_[ed.j]; }
        }
    }
    const double thresh = -p_.eps * (nc ? e / static_cast<double>(nc) : 0.0);
    double gmean = 0.0;
    for (std::size_t i = 0; i < N; ++i) gmean += gsum_[i];
    gmean /= static_cast<double>(N);

    if (!rwin_built_) {
        build_window_table_for(p_.rescue_w, rwin_lo_, rwin_hi_);
        rwin_built_ = true;
    }
    // one rescue exchange per particle per pass (the edge margin budgets exactly one)
    static thread_local std::vector<std::uint8_t> used;
    used.assign(N, 0u);
    if (swap_count_.size() != N) swap_count_.assign(N, 0);

    const std::uint64_t n_cand =
        static_cast<std::uint64_t>(p_.rescue * static_cast<double>(N) + 0.5);
    for (std::uint64_t c = 0; c < n_cand; ++c) {
        const std::uint64_t u = next_u64();
        const std::uint32_t i = static_cast<std::uint32_t>(u % N);
        ++st.rescue_candidates;
        if (used[i]) continue;
        const long ri = static_cast<long>(rank_of_[i]);
        // directed like the main mode-5 channel: over-pressured shrinks, else grows
        const double gi = gsum_[i];
        const bool shrink = (gi < 0.0) ? ((u >> 63) != 0) : (gi > gmean);
        const long a = shrink ? static_cast<long>(rwin_lo_[ri]) : ri + 1;
        const long b = shrink ? ri - 1 : static_cast<long>(rwin_hi_[ri]);
        if (b < a) continue;
        const long r = a + static_cast<long>((u >> 20) % static_cast<std::uint64_t>(b - a + 1));
        const std::uint32_t j = at_rank_[static_cast<std::size_t>(r)];
        if (i == j || used[j]) continue;
        // proxy verdict on this candidate (for the disagreement cross-table)
        const bool i_grows = rad_[i] < rad_[j];
        const double dR = std::fabs(rad_[i] - rad_[j]);
        const double g_grow = i_grows ? gsum_[i] : gsum_[j];
        const double g_shrink = i_grows ? gsum_[j] : gsum_[i];
        const double z2 = static_cast<double>(zc_[i] + zc_[j]);
        const bool proxy_acc = (g_grow - g_shrink) * dR + 0.5 * z2 * dR * dR < thresh;
        const double dE = delta_E(i, j);
        const bool exact_acc = dE < thresh;
        if (exact_acc && !proxy_acc) ++st.rescue_false_rej;
        if (!exact_acc && proxy_acc) ++st.rescue_false_acc;
        if (!exact_acc) continue;
        // apply (same bookkeeping as the exact path)
        std::swap(rad_[i], rad_[j]);
        std::swap(D[i], D[j]);
        std::swap(D0[i], D0[j]);
        const std::uint32_t qi = rank_of_[i], qj = rank_of_[j];
        rank_of_[i] = qj; rank_of_[j] = qi;
        at_rank_[qi] = j; at_rank_[qj] = i;
        touched.push_back(i);
        touched.push_back(j);
        used[i] = used[j] = 1u;
        ++swap_count_[i]; ++swap_count_[j];
        if (p_.event_log_stride > 0 && !event_log_full_ &&
            event_buf_.size() + event_width() <= static_cast<std::size_t>(p_.event_log_cap) * event_width()) {
            event_buf_.push_back(static_cast<double>(current_step_));
            event_buf_.push_back(static_cast<double>(i));
            event_buf_.push_back(static_cast<double>(j));
            event_buf_.push_back((rad_[i] > rad_[j] ? rad_[i] / rad_[j]
                                                    : rad_[j] / rad_[i]) - 1.0);
            event_buf_.push_back(std::fabs(rad_[i] - rad_[j]));
            event_buf_.push_back(-dE);                 // exact relief
            event_buf_.push_back(static_cast<double>(zc_[i]));
            event_buf_.push_back(static_cast<double>(zc_[j]));
            for (std::size_t q = 8; q < event_width(); ++q)
                event_buf_.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        ++st.rescue_accepted;
        st.dE_sum += dE;
        for (std::uint32_t a2 : {i, j}) {
            double g = 0.0; std::uint32_t zn = 0;
            for (std::uint32_t q = adj_off_[a2]; q < adj_off_[a2 + 1]; ++q) {
                const double o = rad_[a2] + rad_[adj_idx_[q]] - adj_r_[q];
                if (o > 0.0) { g += o; ++zn; }
            }
            gsum_[a2] = g; zc_[a2] = zn;
            // neighbours' gsum are stale now, but the rescue pass reads gsum only
            // for direction/proxy bookkeeping, and each particle exchanges at most
            // once per pass -- the exact dE below always decides from adj lists.
        }
    }
    // ---- 3-cycles: trio rotations, all three mutually within rescue_w ----------
    const std::uint64_t n_cyc =
        static_cast<std::uint64_t>(p_.cycle3 * static_cast<double>(N) + 0.5);
    for (std::uint64_t c = 0; c < n_cyc; ++c) {
        const std::uint64_t u = next_u64();
        const std::uint32_t i = static_cast<std::uint32_t>(u % N);
        ++st.cycle3_candidates;
        if (used[i]) continue;
        const long ri = static_cast<long>(rank_of_[i]);
        const long ilo = static_cast<long>(rwin_lo_[ri]), ihi = static_cast<long>(rwin_hi_[ri]);
        if (ihi - ilo < 2) continue;
        const long rj = ilo + static_cast<long>((u >> 20) % static_cast<std::uint64_t>(ihi - ilo + 1));
        if (rj == ri) continue;
        const std::uint32_t j = at_rank_[static_cast<std::size_t>(rj)];
        if (used[j]) continue;
        // k mutually within rescue_w of BOTH i and j: intersect the rank windows
        const long klo = std::max(ilo, static_cast<long>(rwin_lo_[rj]));
        const long khi = std::min(ihi, static_cast<long>(rwin_hi_[rj]));
        if (khi < klo) continue;
        const long rk = klo + static_cast<long>((u >> 41) % static_cast<std::uint64_t>(khi - klo + 1));
        if (rk == ri || rk == rj) continue;
        const std::uint32_t k = at_rank_[static_cast<std::size_t>(rk)];
        if (used[k]) continue;
        // try both orientations; take the first that descends
        std::uint32_t a = i, b = j, cc = k;
        double dE = delta_E3(a, b, cc);
        if (!(dE < thresh)) {
            std::swap(b, cc);
            dE = delta_E3(a, b, cc);
            if (!(dE < thresh)) continue;
        }
        // pairwise-blocked instrumentation: were all three pair swaps uphill?
        const bool blocked = !(delta_E(i, j) < thresh) && !(delta_E(j, k) < thresh) &&
                             !(delta_E(i, k) < thresh);
        // apply rotation a<-D(b), b<-D(cc), cc<-D(a)
        const double Da = D[a], D0a = D0[a], Ra = rad_[a];
        D[a] = D[b];   D0[a] = D0[b];   rad_[a] = rad_[b];
        D[b] = D[cc];  D0[b] = D0[cc];  rad_[b] = rad_[cc];
        D[cc] = Da;    D0[cc] = D0a;    rad_[cc] = Ra;
        const std::uint32_t qa = rank_of_[a], qb = rank_of_[b], qc = rank_of_[cc];
        // a now holds b's D0 (rank qb), b holds cc's (rank qc), cc holds a's (rank qa)
        rank_of_[a] = qb; rank_of_[b] = qc; rank_of_[cc] = qa;
        at_rank_[qb] = a; at_rank_[qc] = b; at_rank_[qa] = cc;
        touched.push_back(a); touched.push_back(b); touched.push_back(cc);
        used[a] = used[b] = used[cc] = 1u;
        ++swap_count_[a]; ++swap_count_[b]; ++swap_count_[cc];
        ++st.cycle3_accepted;
        if (blocked) ++st.cycle3_blocked;
        st.dE_sum += dE;
        for (std::uint32_t m : {a, b, cc}) {
            double g = 0.0; std::uint32_t zn = 0;
            for (std::uint32_t q = adj_off_[m]; q < adj_off_[m + 1]; ++q) {
                const double o = rad_[m] + rad_[adj_idx_[q]] - adj_r_[q];
                if (o > 0.0) { g += o; ++zn; }
            }
            gsum_[m] = g; zc_[m] = zn;
        }
    }
    rescue_cand_total_ += st.rescue_candidates;
    rescue_acc_total_ += st.rescue_accepted;
    rescue_frej_total_ += st.rescue_false_rej;
    rescue_facc_total_ += st.rescue_false_acc;
    cycle3_cand_total_ += st.cycle3_candidates;
    cycle3_acc_total_ += st.cycle3_accepted;
    cycle3_blocked_total_ += st.cycle3_blocked;
}

}  // namespace rcpgenerator
