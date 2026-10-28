/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_SPECTRAL_TRACKER_H
#define HEIPROMAP_SPECTRAL_TRACKER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "aligned_array.h"
#include "mapping.h"
#include "profiler.h"

namespace HeiProMap {

    // Effort of the re-smoothing after each prolongation: rounds of inexact inverse iteration, each with `depth` CG
    // iterations. Levels with more than big_n vertices get their own, cheaper settings: they dominate the run time,
    // never reach the residual tolerance anyway, and a single round gave the better partition on the 10M-row test.
    struct SpectralSmoothing {
        int rounds = 3;
        int depth = 24;
        int big_rounds = 1;
        int big_depth = 8;
        vertex_t big_n = 500000;
        int threads = 1; // threads of the vector operations and Laplacian products; the result does not depend on it
    };

    // Bottom Laplacian subspace of the current graph level, carried along the
    // solver's own coarsening hierarchy (multilevel inverse iteration as in
    // Barnard & Simon 1994), plus recursive coordinate bisection onto the grid.
    class SpectralTracker {
    public:
        // A vector of the current number of vertices (n_) whose elements are not initialized: every use writes it
        // before reading it. Its size() is rounded up, so the loops use n_.
        using Vec = AlignedArray<f64>;

        static constexpr size_t kSubspace = 3;
        static constexpr vertex_t kDenseMaxN = 4096;
        static constexpr int kCoarseRounds = 20;
        static constexpr int kFallbackRounds = 30;
        static constexpr int kCgDepth = 24; // work budget of the CG fallback for the coarsest level
        static constexpr f64 kResidualTol = 0.03;

        void set_smoothing(const SpectralSmoothing &smoothing) { smoothing_ = smoothing; }

        bool initialized() const { return n_ > 0; }

        // The kSubspace mean-free orthonormal vectors, X[0] smoothest.
        const std::vector<Vec> &basis() const { return X; }

        vertex_t size() const { return static_cast<vertex_t>(n_); }

        void init(const CSRGraph &g, u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "init");
            n_ = g.n;
            X = make_vecs(kSubspace, n_);
            std::mt19937_64 rng(seed ^ 0x5BEC7A1ull);
            std::uniform_real_distribution<f64> unif(-0.5, 0.5);
            for (auto &x: X) for (size_t i = 0; i < n_; ++i) x[i] = unif(rng);

            if (g.n <= kDenseMaxN) {
                dense_solve(g);
            } else {
                resmooth(g, kFallbackRounds, 4 * kCgDepth);
            }
        }

        // Piecewise-constant prolongation onto the finer graph, then re-smoothing.
        void prolong(const Mapping &m, const CSRGraph &fine) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "prolong");
            const size_t n = fine.n;
            std::vector<Vec> F = make_vecs(kSubspace, n);
            for (size_t k = 0; k < kSubspace; ++k) {
                #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
                for (size_t u = 0; u < n; ++u) {
                    F[k][u] = X[k][m.get(static_cast<vertex_t>(u))];
                }
            }
            X.swap(F);
            n_ = n;
            if (fine.n > smoothing_.big_n) {
                resmooth(fine, smoothing_.big_rounds, smoothing_.big_depth);
            } else {
                resmooth(fine, smoothing_.rounds, smoothing_.depth);
            }
        }

        // Recursive coordinate bisection onto the W x H grid described by topo
        // (cell id = row * W + col, capacity = topo vertex weight). Column cuts use
        // the smoothest direction, row cuts the second, as in SpectralEmbeddingMapper.
        void place(const CSRGraph &g, const CSRGraph &topo, size_t W, std::vector<partition_t> &out) const {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "place");
            const size_t k = topo.n;
            const size_t H = (k + W - 1) / W;
            const Vec &cx = W >= H ? X[0] : X[1];
            const Vec &cy = W >= H ? X[1] : X[0];

            // 2D prefix sums of cell capacities for O(1) rectangle sums.
            std::vector<f64> cap_ps((W + 1) * (H + 1), 0.0);
            for (size_t r = 0; r < H; ++r) {
                for (size_t c = 0; c < W; ++c) {
                    const size_t id = r * W + c;
                    const f64 w = id < k ? static_cast<f64>(topo.v_weight(static_cast<vertex_t>(id))) : 0.0;
                    cap_ps[(r + 1) * (W + 1) + (c + 1)] = w + cap_ps[r * (W + 1) + (c + 1)]
                                                          + cap_ps[(r + 1) * (W + 1) + c] - cap_ps[r * (W + 1) + c];
                }
            }

            // A rectangle of cells and the vertices perm[lo, hi) that are placed in it.
            struct Node {
                size_t lo, hi, col0, row0, w, h;
            };

            out.assign(g.n, 0);
            std::vector<vertex_t> perm(g.n);
            std::iota(perm.begin(), perm.end(), static_cast<vertex_t>(0));

            std::vector<Node> stack{{0, static_cast<size_t>(g.n), 0, 0, W, H}};
            while (!stack.empty()) {
                const Node nd = stack.back();
                stack.pop_back();

                if (nd.w * nd.h == 1) {
                    const auto cell = static_cast<partition_t>(nd.row0 * W + nd.col0);
                    for (size_t i = nd.lo; i < nd.hi; ++i) out[perm[i]] = cell;
                    continue;
                }

                const bool split_cols = nd.w >= nd.h;
                Node a{}, b{};
                if (split_cols) {
                    const size_t wA = nd.w / 2;
                    a = {0, 0, nd.col0, nd.row0, wA, nd.h};
                    b = {0, 0, nd.col0 + wA, nd.row0, nd.w - wA, nd.h};
                } else {
                    const size_t hA = nd.h / 2;
                    a = {0, 0, nd.col0, nd.row0, nd.w, hA};
                    b = {0, 0, nd.col0, nd.row0 + hA, nd.w, nd.h - hA};
                }

                const KeyLess less{split_cols ? &cx : &cy};
                std::sort(perm.begin() + nd.lo, perm.begin() + nd.hi, less);

                const f64 capA = rect_capacity(cap_ps, W, a.col0, a.row0, a.w, a.h);
                const f64 cap = capA + rect_capacity(cap_ps, W, b.col0, b.row0, b.w, b.h);
                f64 total_w = 0.0;
                for (size_t i = nd.lo; i < nd.hi; ++i) total_w += static_cast<f64>(g.v_weight(perm[i]));
                const f64 target = cap > 0.0 ? total_w * capA / cap : 0.0;

                // Split where the prefix weight is closest to side A's capacity share.
                size_t split = nd.lo;
                f64 prefix = 0.0, best_gap = std::abs(target);
                for (size_t i = nd.lo; i < nd.hi; ++i) {
                    prefix += static_cast<f64>(g.v_weight(perm[i]));
                    const f64 gap = std::abs(prefix - target);
                    if (gap < best_gap) {
                        best_gap = gap;
                        split = i + 1;
                    }
                }
                const size_t m = nd.hi - nd.lo;
                if (m >= 2) {
                    if (split == nd.lo && capA > 0.0) split = nd.lo + 1;
                    if (split == nd.hi && cap > capA) split = nd.hi - 1;
                }

                a.lo = nd.lo;
                a.hi = split;
                b.lo = split;
                b.hi = nd.hi;
                stack.push_back(b);
                stack.push_back(a);
            }
        }

        // Same objective as get_message_cost(), evaluated from a plain partition array.
        template<typename PartT, typename OracleT>
        static f64 message_cost(const CSRGraph &g, const PartT &part, const OracleT &oracle,
                                f64 alpha_end, f64 beta_end, f64 alpha_link, f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "message_cost");
            const u64 k = oracle.current_k();
            std::vector<std::pair<u64, f64> > pairs;
            pairs.reserve(g.m);
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t bu = part[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t bv = part[g.edges_v[j]];
                    if (bu != bv) pairs.emplace_back(static_cast<u64>(bu) * k + bv, static_cast<f64>(g.edge_w(j)));
                }
            }
            std::sort(pairs.begin(), pairs.end());

            f64 cost = 0.0;
            for (size_t i = 0; i < pairs.size();) {
                const u64 key = pairs[i].first;
                f64 S = 0.0;
                for (; i < pairs.size() && pairs[i].first == key; ++i) S += pairs[i].second;
                const f64 h = static_cast<f64>(oracle.get(static_cast<partition_t>(key / k), static_cast<partition_t>(key % k)));
                cost += alpha_end + beta_end * S + alpha_link * h + beta_link * h * S;
            }
            return cost;
        }

    private:
        std::vector<Vec> X; // kSubspace mean-free orthonormal vectors, X[0] smoothest
        size_t n_ = 0;      // their length: the number of vertices of the current level
        SpectralSmoothing smoothing_;

        // Reductions add the sums of chunks of fixed size in order, whatever the number of threads, so every result
        // is the same for any number of threads. The other loops write independent elements.
        static constexpr size_t kChunk = size_t(1) << 14;
        static constexpr size_t kParallelMin = size_t(1) << 16; // loops over fewer elements run on one thread

        // Threads of a loop over n elements.
        int threads_for(const size_t n) const { return n >= kParallelMin ? std::max(1, smoothing_.threads) : 1; }

        static size_t num_chunks(const size_t n) { return (n + kChunk - 1) / kChunk; }

        // Orders the vertices by one coordinate, ties by id.
        struct KeyLess {
            const Vec *key;

            bool operator()(const vertex_t x, const vertex_t y) const {
                return (*key)[x] != (*key)[y] ? (*key)[x] < (*key)[y] : x < y;
            }
        };

        // Capacity of the w x h cells with the top left cell (c0, r0), from the 2D prefix sums.
        static f64 rect_capacity(const std::vector<f64> &ps, const size_t W,
                                 const size_t c0, const size_t r0, const size_t w, const size_t h) {
            return ps[(r0 + h) * (W + 1) + (c0 + w)] - ps[r0 * (W + 1) + (c0 + w)]
                   - ps[(r0 + h) * (W + 1) + c0] + ps[r0 * (W + 1) + c0];
        }

        // count vectors of n elements each, not initialized.
        static std::vector<Vec> make_vecs(const size_t count, const size_t n) {
            std::vector<Vec> v(count);
            for (auto &e: v) { e.initialize(n); }
            return v;
        }

        // Weighted degree of every vertex.
        Vec degrees(const CSRGraph &g) const {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "degrees");
            Vec deg;
            deg.initialize(g.n);
            #pragma omp parallel for num_threads(threads_for(g.n)) schedule(static)
            for (size_t u = 0; u < g.n; ++u) {
                f64 d = 0.0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) d += static_cast<f64>(g.edge_w(j));
                deg[u] = d;
            }
            return deg;
        }

        // Compact read-only copy of a graph for the bandwidth-bound Laplacian products: 32-bit neighbor ids and,
        // if the graph is not unweighted, 32-bit float weights (kept as doubles if some weight is not exactly
        // representable). Weights are converted back to f64 exactly, so products equal those on the original graph.
        struct CompactCsr {
            vertex_t n = 0;
            const size_t *off = nullptr;
            AlignedArray<u32> adj;
            AlignedArray<f32> wf;
            AlignedArray<f64> wd;
            bool uniform = true;
            bool weights_f32 = false; // otherwise wd holds the weights (if the graph is not unweighted)

            CompactCsr(const CSRGraph &g, const int threads) {
                HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "compact_csr");
                n = g.n;
                off = &g.neighborhoods[0];
                if (g.n >= (static_cast<vertex_t>(1) << 32)) {
                    std::cout << "Error: spectral solver supports at most 2^32 - 1 vertices" << std::endl;
                    exit(EXIT_FAILURE);
                }
                const size_t m = g.neighborhoods[g.n];
                const int th = m >= kParallelMin ? std::max(1, threads) : 1;
                adj.initialize(m);
                #pragma omp parallel for num_threads(th) schedule(static)
                for (size_t j = 0; j < m; ++j) adj[j] = static_cast<u32>(g.edges_v[j]);

                uniform = g.uniform_e_weights;
                if (uniform) return;

                constexpr weight_t kExactF32 = static_cast<weight_t>(1) << 24;
                size_t n_inexact = 0;
                #pragma omp parallel for num_threads(th) schedule(static) reduction(+:n_inexact)
                for (size_t j = 0; j < m; ++j) n_inexact += (g.edges_w[j] < 0 || g.edges_w[j] > kExactF32);
                weights_f32 = n_inexact == 0;
                if (weights_f32) {
                    wf.initialize(m);
                    #pragma omp parallel for num_threads(th) schedule(static)
                    for (size_t j = 0; j < m; ++j) wf[j] = static_cast<f32>(g.edges_w[j]);
                } else {
                    wd.initialize(m);
                    #pragma omp parallel for num_threads(th) schedule(static)
                    for (size_t j = 0; j < m; ++j) wd[j] = static_cast<f64>(g.edges_w[j]);
                }
            }

            f64 weight(const size_t j) const {
                if (uniform) return 1.0;
                return weights_f32 ? static_cast<f64>(wf[j]) : wd[j];
            }
        };

        // out = L in for one vector.
        void apply_L(const CompactCsr &c, const Vec &deg, const Vec &in, Vec &out) const {
            #pragma omp parallel for num_threads(threads_for(c.n)) schedule(dynamic, 2048)
            for (size_t u = 0; u < c.n; ++u) {
                f64 s = deg[u] * in[u];
                for (size_t j = c.off[u]; j < c.off[u + 1]; ++j) s -= c.weight(j) * in[c.adj[j]];
                out[u] = s;
            }
        }

        // out_k = L in_k for the kSubspace vectors at once: the matrix is streamed once instead of kSubspace times.
        // Each output is computed with the same operations in the same order as apply_L.
        void apply_L_block(const CompactCsr &c, const Vec &deg,
                           const std::array<const Vec *, kSubspace> &in,
                           const std::array<Vec *, kSubspace> &out) const {
            static_assert(kSubspace == 3, "apply_L_block is written for three vectors");
            const f64 *x0 = in[0]->get_ptr(), *x1 = in[1]->get_ptr(), *x2 = in[2]->get_ptr();
            f64 *y0 = out[0]->get_ptr(), *y1 = out[1]->get_ptr(), *y2 = out[2]->get_ptr();
            #pragma omp parallel for num_threads(threads_for(c.n)) schedule(dynamic, 2048)
            for (size_t u = 0; u < c.n; ++u) {
                f64 s0 = deg[u] * x0[u], s1 = deg[u] * x1[u], s2 = deg[u] * x2[u];
                for (size_t j = c.off[u]; j < c.off[u + 1]; ++j) {
                    const f64 wj = c.weight(j);
                    const u32 v = c.adj[j];
                    s0 -= wj * x0[v];
                    s1 -= wj * x1[v];
                    s2 -= wj * x2[v];
                }
                y0[u] = s0;
                y1[u] = s1;
                y2[u] = s2;
            }
        }

        f64 dot(const Vec &a, const Vec &b) const {
            const size_t n = n_;
            const size_t chunks = num_chunks(n);
            std::vector<f64> partial(chunks);
            #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
            for (size_t ch = 0; ch < chunks; ++ch) {
                const size_t lo = ch * kChunk, hi = std::min(n, lo + kChunk);
                f64 s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                size_t i = lo;
                for (; i + 4 <= hi; i += 4) {
                    s0 += a[i] * b[i];
                    s1 += a[i + 1] * b[i + 1];
                    s2 += a[i + 2] * b[i + 2];
                    s3 += a[i + 3] * b[i + 3];
                }
                for (; i < hi; ++i) s0 += a[i] * b[i];
                partial[ch] = (s0 + s1) + (s2 + s3);
            }
            f64 s = 0.0;
            for (const f64 p: partial) s += p;
            return s;
        }

        void subtract_mean(Vec &v) const {
            if (n_ == 0) return;
            const size_t n = n_;
            const size_t chunks = num_chunks(n);
            std::vector<f64> partial(chunks);
            #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
            for (size_t ch = 0; ch < chunks; ++ch) {
                const size_t lo = ch * kChunk, hi = std::min(n, lo + kChunk);
                f64 s = 0.0;
                for (size_t i = lo; i < hi; ++i) s += v[i];
                partial[ch] = s;
            }
            f64 sum = 0.0;
            for (const f64 p: partial) sum += p;
            const f64 mean = sum / static_cast<f64>(n);
            #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
            for (size_t i = 0; i < n; ++i) v[i] -= mean;
        }

        // y += a * x
        void axpy(Vec &y, const f64 a, const Vec &x) const {
            const size_t n = n_;
            #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
            for (size_t i = 0; i < n; ++i) y[i] += a * x[i];
        }

        // y *= a
        void scale(Vec &y, const f64 a) const {
            const size_t n = n_;
            #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
            for (size_t i = 0; i < n; ++i) y[i] *= a;
        }

        // Removes the constant vector and the components along X[0..k) from X[k].
        void project_out(const size_t k) {
            subtract_mean(X[k]);
            for (size_t j = 0; j < k; ++j) axpy(X[k], -dot(X[k], X[j]), X[j]);
        }

        // Modified Gram-Schmidt against the constant vector and each other;
        // a collapsed vector is re-seeded.
        void orthonormalise(u64 salt) {
            for (size_t k = 0; k < X.size(); ++k) {
                project_out(k);
                f64 nk = std::sqrt(dot(X[k], X[k]));
                if (!(nk > 1e-30)) {
                    std::mt19937_64 rng(salt * 31 + k);
                    std::uniform_real_distribution<f64> unif(-0.5, 0.5);
                    for (size_t i = 0; i < n_; ++i) X[k][i] = unif(rng);
                    project_out(k);
                    nk = std::sqrt(dot(X[k], X[k]));
                    if (!(nk > 0.0)) nk = 1.0;
                }
                scale(X[k], 1.0 / nk);
            }
        }

        // Cyclic Jacobi on a tiny symmetric s x s matrix; V[i*s+j] = component i of eigenvector j.
        static void small_sym_eig(size_t s, std::vector<f64> &B, std::vector<f64> &V) {
            V.assign(s * s, 0.0);
            for (size_t i = 0; i < s; ++i) V[i * s + i] = 1.0;
            for (int sweep = 0; sweep < 32; ++sweep) {
                f64 off = 0.0;
                for (size_t p = 0; p < s; ++p)
                    for (size_t q = p + 1; q < s; ++q) off += B[p * s + q] * B[p * s + q];
                if (off <= 1e-30) break;
                for (size_t p = 0; p < s; ++p) {
                    for (size_t q = p + 1; q < s; ++q) {
                        const f64 bpq = B[p * s + q];
                        if (bpq == 0.0) continue;
                        const f64 theta = (B[q * s + q] - B[p * s + p]) / (2.0 * bpq);
                        const f64 t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                        const f64 c = 1.0 / std::sqrt(t * t + 1.0), sn = t * c;
                        for (size_t i = 0; i < s; ++i) {
                            const f64 bip = B[i * s + p], biq = B[i * s + q];
                            B[i * s + p] = c * bip - sn * biq;
                            B[i * s + q] = sn * bip + c * biq;
                        }
                        for (size_t i = 0; i < s; ++i) {
                            const f64 bpi = B[p * s + i], bqi = B[q * s + i];
                            B[p * s + i] = c * bpi - sn * bqi;
                            B[q * s + i] = sn * bpi + c * bqi;
                            const f64 vip = V[i * s + p], viq = V[i * s + q];
                            V[i * s + p] = c * vip - sn * viq;
                            V[i * s + q] = sn * vip + c * viq;
                        }
                    }
                }
            }
        }

        // Rotate the (orthonormal) basis into Ritz vectors ordered by ascending eigenvalue.
        void rayleigh_ritz(const CompactCsr &g, const Vec &deg) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "rayleigh_ritz");
            const size_t s = X.size();
            std::vector<Vec> LX = make_vecs(s, g.n);
            apply_L_block(g, deg, {&X[0], &X[1], &X[2]}, {&LX[0], &LX[1], &LX[2]});

            std::vector<f64> B(s * s);
            for (size_t i = 0; i < s; ++i)
                for (size_t j = i; j < s; ++j) B[i * s + j] = B[j * s + i] = dot(X[i], LX[j]);

            std::vector<f64> V;
            small_sym_eig(s, B, V);

            // Eigenvectors by ascending eigenvalue (ties by index): insertion sort of s indices.
            std::vector<size_t> ord(s);
            std::iota(ord.begin(), ord.end(), static_cast<size_t>(0));
            for (size_t i = 1; i < s; ++i) {
                for (size_t j = i; j > 0; --j) {
                    const size_t a = ord[j], b = ord[j - 1];
                    const f64 ea = B[a * s + a], eb = B[b * s + b];
                    if (!(ea < eb || (ea == eb && a < b))) break;
                    std::swap(ord[j], ord[j - 1]);
                }
            }

            std::vector<Vec> NX = make_vecs(s, g.n);
            #pragma omp parallel for num_threads(threads_for(g.n)) schedule(static)
            for (size_t r = 0; r < g.n; ++r) {
                for (size_t k = 0; k < s; ++k) {
                    f64 acc = 0.0;
                    for (size_t j = 0; j < s; ++j) acc += X[j][r] * V[j * s + ord[k]];
                    NX[k][r] = acc;
                }
            }
            X.swap(NX);
        }

        // Fixed-depth CG on L x_k = b_k for the kSubspace right-hand sides b_k = x_k (b_k mean-free), run in
        // lockstep so that all vectors share one matrix pass per iteration. Every vector keeps its own step sizes
        // and stops on its own; the depth is a work budget, not a tolerance. On return x_k holds the solution.
        void cg_solve_block(const CompactCsr &g, const Vec &deg, int depth) {
            const size_t n = g.n;
            constexpr size_t S = kSubspace;
            std::array<Vec, S> r, p, Ap;
            std::array<f64, S> rsq{};
            std::array<bool, S> active{};
            for (size_t k = 0; k < S; ++k) {
                r[k].initialize(n);
                p[k].initialize(n);
                Ap[k].initialize(n);
                #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
                for (size_t i = 0; i < n; ++i) {
                    r[k][i] = X[k][i];
                    p[k][i] = X[k][i];
                    X[k][i] = 0.0;
                }
                rsq[k] = dot(r[k], r[k]);
                active[k] = rsq[k] > 0.0;
            }
            for (int it = 0; it < depth; ++it) {
                if (!(active[0] || active[1] || active[2])) break;
                apply_L_block(g, deg, {&p[0], &p[1], &p[2]}, {&Ap[0], &Ap[1], &Ap[2]});
                for (size_t k = 0; k < S; ++k) {
                    if (!active[k]) continue;
                    const f64 pAp = dot(p[k], Ap[k]);
                    if (!(pAp > 0.0)) {
                        active[k] = false;
                        continue;
                    }
                    const f64 alpha = rsq[k] / pAp;
                    #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
                    for (size_t i = 0; i < n; ++i) {
                        X[k][i] += alpha * p[k][i];
                        r[k][i] -= alpha * Ap[k][i];
                    }
                    const f64 nrsq = dot(r[k], r[k]);
                    const f64 beta = nrsq / rsq[k];
                    rsq[k] = nrsq;
                    #pragma omp parallel for num_threads(threads_for(n)) schedule(static)
                    for (size_t i = 0; i < n; ++i) p[k][i] = r[k][i] + beta * p[k][i];
                }
            }
        }

        // Relative residual ||L x - lambda x|| / lambda of the leading Ritz pair.
        f64 leading_residual(const CompactCsr &g, const Vec &deg, Vec &scratch) const {
            apply_L(g, deg, X[0], scratch);
            const f64 lam = dot(X[0], scratch);
            const f64 r2 = dot(scratch, scratch) - lam * lam;
            return std::sqrt(std::max(0.0, r2)) / std::max(lam, 1e-12);
        }

        // Rounds of inexact inverse iteration + Rayleigh-Ritz; returns {rounds, residual}.
        std::pair<int, f64> resmooth(const CSRGraph &g, int budget, int cg_depth) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "resmooth");
            const Vec deg = degrees(g);
            const CompactCsr compact(g, smoothing_.threads);
            Vec scratch;
            scratch.initialize(g.n);
            orthonormalise(0xA110u);
            int rounds = 0;
            f64 res = 0.0;
            for (int r = 0; r < budget; ++r) {
                {
                    HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "cg");
                    cg_solve_block(compact, deg, cg_depth);
                }
                {
                    HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "orthonormalise");
                    for (auto &xk: X) subtract_mean(xk);
                    orthonormalise(static_cast<u64>(r));
                }
                rayleigh_ritz(compact, deg);
                {
                    HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "orthonormalise");
                    orthonormalise(static_cast<u64>(r) + 0x8000u);
                }
                {
                    HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "residual");
                    res = leading_residual(compact, deg, scratch);
                }
                ++rounds;
                if (res <= kResidualTol) break;
            }
            return {rounds, res};
        }

        // Solves (L + delta*I) x = b in place with the Cholesky factor f (lower triangle, row-major, n x n).
        static void cholesky_solve(const std::vector<f64> &f, const size_t n, Vec &x) {
            for (size_t i = 0; i < n; ++i) {
                f64 s = x[i];
                for (size_t k = 0; k < i; ++k) s -= f[i * n + k] * x[k];
                x[i] = s / f[i * n + i];
            }
            for (size_t i = n; i-- > 0;) {
                f64 s = x[i];
                for (size_t k = i + 1; k < n; ++k) s -= f[k * n + i] * x[k];
                x[i] = s / f[i * n + i];
            }
        }

        // Exact inverse subspace iteration with a dense Cholesky factor of L + delta*I.
        void dense_solve(const CSRGraph &g) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "dense_solve");
            const size_t n = g.n;
            const Vec deg = degrees(g);
            const f64 maxdeg = n == 0 ? 0.0 : *std::max_element(deg.begin(), deg.begin() + n);
            const f64 delta = std::max(1e-300, 1e-8 * maxdeg);

            std::vector<f64> f(n * n, 0.0);
            for (vertex_t v = 0; v < g.n; ++v) {
                f[v * n + v] += deg[v] + delta;
                for (size_t j = g.neighborhoods[v]; j < g.neighborhoods[v + 1]; ++j) {
                    const vertex_t u = g.edges_v[j];
                    if (u < v) f[v * n + u] -= static_cast<f64>(g.edge_w(j));
                }
            }
            for (size_t j = 0; j < n; ++j) {
                f64 d = f[j * n + j];
                for (size_t k = 0; k < j; ++k) d -= f[j * n + k] * f[j * n + k];
                if (!(d > 0.0)) {
                    std::cout << "  [spectral] dense Cholesky failed, falling back to CG" << std::endl;
                    resmooth(g, kFallbackRounds, 4 * kCgDepth);
                    return;
                }
                const f64 dj = std::sqrt(d);
                f[j * n + j] = dj;
                for (size_t i = j + 1; i < n; ++i) {
                    f64 s = f[i * n + j];
                    for (size_t k = 0; k < j; ++k) s -= f[i * n + k] * f[j * n + k];
                    f[i * n + j] = s / dj;
                }
            }

            for (int r = 0; r < kCoarseRounds; ++r) {
                for (auto &xk: X) {
                    cholesky_solve(f, n, xk);
                    subtract_mean(xk);
                }
                orthonormalise(static_cast<u64>(r));
            }
            rayleigh_ritz(CompactCsr(g, smoothing_.threads), deg);
            orthonormalise(0x5EEDu);
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_SPECTRAL_TRACKER_H
