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
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "mapping.h"
#include "profiler.h"

namespace HeiProMap {

    // Bottom Laplacian subspace of the current graph level, carried along the
    // solver's own coarsening hierarchy (multilevel inverse iteration as in
    // Barnard & Simon 1994), plus recursive coordinate bisection onto the grid.
    class SpectralTracker {
    public:
        static constexpr size_t kSubspace = 3;
        static constexpr vertex_t kDenseMaxN = 4096;
        static constexpr int kCoarseRounds = 20;
        static constexpr int kFallbackRounds = 30;
        static constexpr int kMaxRounds = 3;
        static constexpr int kCgDepth = 24;
        static constexpr f64 kResidualTol = 0.03;

        bool initialized() const { return !X.empty() && !X[0].empty(); }

        vertex_t size() const { return initialized() ? static_cast<vertex_t>(X[0].size()) : 0; }

        void init(const CSRGraph &g, u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "init");
            const size_t n = g.n;
            X.assign(kSubspace, std::vector<f64>(n));
            std::mt19937_64 rng(seed ^ 0x5BEC7A1ull);
            std::uniform_real_distribution<f64> unif(-0.5, 0.5);
            for (auto &x: X) for (auto &e: x) e = unif(rng);

            if (g.n <= kDenseMaxN) {
                dense_solve(g);
            } else {
                resmooth(g, kFallbackRounds, 4 * kCgDepth);
            }
        }

        // Piecewise-constant prolongation onto the finer graph, then re-smoothing.
        void prolong(const Mapping &m, const CSRGraph &fine) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "prolong");
            std::vector<std::vector<f64>> F(kSubspace, std::vector<f64>(fine.n));
            for (size_t k = 0; k < kSubspace; ++k) {
                for (vertex_t u = 0; u < fine.n; ++u) F[k][u] = X[k][m.get(u)];
            }
            X.swap(F);
            resmooth(fine, kMaxRounds, kCgDepth);
        }

        // Recursive coordinate bisection onto the W x H grid described by topo
        // (cell id = row * W + col, capacity = topo vertex weight). Column cuts use
        // the smoothest direction, row cuts the second, as in SpectralEmbeddingMapper.
        void place(const CSRGraph &g, const CSRGraph &topo, size_t W, std::vector<partition_t> &out) const {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "place");
            const size_t k = topo.n;
            const size_t H = (k + W - 1) / W;
            const std::vector<f64> &cx = W >= H ? X[0] : X[1];
            const std::vector<f64> &cy = W >= H ? X[1] : X[0];

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
            auto rect_cap = [&](size_t c0, size_t r0, size_t w, size_t h) {
                return cap_ps[(r0 + h) * (W + 1) + (c0 + w)] - cap_ps[r0 * (W + 1) + (c0 + w)]
                       - cap_ps[(r0 + h) * (W + 1) + c0] + cap_ps[r0 * (W + 1) + c0];
            };

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

                const std::vector<f64> &key = split_cols ? cx : cy;
                std::sort(perm.begin() + nd.lo, perm.begin() + nd.hi, [&](vertex_t x, vertex_t y) {
                    return key[x] != key[y] ? key[x] < key[y] : x < y;
                });

                const f64 capA = rect_cap(a.col0, a.row0, a.w, a.h);
                const f64 cap = capA + rect_cap(b.col0, b.row0, b.w, b.h);
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
        std::vector<std::vector<f64>> X; // kSubspace mean-free orthonormal vectors, X[0] smoothest

        static std::vector<f64> degrees(const CSRGraph &g) {
            std::vector<f64> deg(g.n, 0.0);
            for (vertex_t u = 0; u < g.n; ++u) {
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) deg[u] += static_cast<f64>(g.edge_w(j));
            }
            return deg;
        }

        static void apply_L(const CSRGraph &g, const std::vector<f64> &deg, const std::vector<f64> &in, std::vector<f64> &out) {
            for (vertex_t u = 0; u < g.n; ++u) {
                f64 s = deg[u] * in[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    s -= static_cast<f64>(g.edge_w(j)) * in[g.edges_v[j]];
                }
                out[u] = s;
            }
        }

        static f64 dot(const std::vector<f64> &a, const std::vector<f64> &b) {
            f64 s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            const size_t n = a.size();
            size_t i = 0;
            for (; i + 4 <= n; i += 4) {
                s0 += a[i] * b[i];
                s1 += a[i + 1] * b[i + 1];
                s2 += a[i + 2] * b[i + 2];
                s3 += a[i + 3] * b[i + 3];
            }
            for (; i < n; ++i) s0 += a[i] * b[i];
            return (s0 + s1) + (s2 + s3);
        }

        static void subtract_mean(std::vector<f64> &v) {
            if (v.empty()) return;
            const f64 mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<f64>(v.size());
            for (auto &e: v) e -= mean;
        }

        // Modified Gram-Schmidt against the constant vector and each other;
        // a collapsed vector is re-seeded.
        void orthonormalise(u64 salt) {
            for (size_t k = 0; k < X.size(); ++k) {
                auto &xk = X[k];
                auto project_out = [&]() {
                    subtract_mean(xk);
                    for (size_t j = 0; j < k; ++j) {
                        const f64 p = dot(xk, X[j]);
                        for (size_t i = 0; i < xk.size(); ++i) xk[i] -= p * X[j][i];
                    }
                };
                project_out();
                f64 nk = std::sqrt(dot(xk, xk));
                if (!(nk > 1e-30)) {
                    std::mt19937_64 rng(salt * 31 + k);
                    std::uniform_real_distribution<f64> unif(-0.5, 0.5);
                    for (auto &e: xk) e = unif(rng);
                    project_out();
                    nk = std::sqrt(dot(xk, xk));
                    if (!(nk > 0.0)) nk = 1.0;
                }
                const f64 inv = 1.0 / nk;
                for (auto &e: xk) e *= inv;
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
        void rayleigh_ritz(const CSRGraph &g, const std::vector<f64> &deg) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "rayleigh_ritz");
            const size_t s = X.size();
            std::vector<std::vector<f64>> LX(s, std::vector<f64>(g.n));
            for (size_t k = 0; k < s; ++k) apply_L(g, deg, X[k], LX[k]);

            std::vector<f64> B(s * s);
            for (size_t i = 0; i < s; ++i)
                for (size_t j = i; j < s; ++j) B[i * s + j] = B[j * s + i] = dot(X[i], LX[j]);

            std::vector<f64> V;
            small_sym_eig(s, B, V);
            std::vector<size_t> ord(s);
            std::iota(ord.begin(), ord.end(), static_cast<size_t>(0));
            std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
                return B[a * s + a] != B[b * s + b] ? B[a * s + a] < B[b * s + b] : a < b;
            });

            std::vector<std::vector<f64>> NX(s, std::vector<f64>(g.n, 0.0));
            for (vertex_t r = 0; r < g.n; ++r) {
                for (size_t k = 0; k < s; ++k) {
                    f64 acc = 0.0;
                    for (size_t j = 0; j < s; ++j) acc += X[j][r] * V[j * s + ord[k]];
                    NX[k][r] = acc;
                }
            }
            X.swap(NX);
        }

        // Fixed-depth CG on L x = b (b mean-free): the depth is a work budget, not a tolerance.
        static void cg_solve(const CSRGraph &g, const std::vector<f64> &deg, const std::vector<f64> &b,
                             std::vector<f64> &x, int depth) {
            const size_t n = b.size();
            std::vector<f64> r(b), p(b), Ap(n);
            std::fill(x.begin(), x.end(), 0.0);
            f64 rsq = dot(r, r);
            if (!(rsq > 0.0)) return;
            for (int it = 0; it < depth; ++it) {
                apply_L(g, deg, p, Ap);
                const f64 pAp = dot(p, Ap);
                if (!(pAp > 0.0)) break;
                const f64 alpha = rsq / pAp;
                for (size_t i = 0; i < n; ++i) {
                    x[i] += alpha * p[i];
                    r[i] -= alpha * Ap[i];
                }
                const f64 nrsq = dot(r, r);
                const f64 beta = nrsq / rsq;
                rsq = nrsq;
                for (size_t i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
            }
        }

        // Relative residual ||L x - lambda x|| / lambda of the leading Ritz pair.
        f64 leading_residual(const CSRGraph &g, const std::vector<f64> &deg, std::vector<f64> &scratch) const {
            apply_L(g, deg, X[0], scratch);
            const f64 lam = dot(X[0], scratch);
            const f64 r2 = dot(scratch, scratch) - lam * lam;
            return std::sqrt(std::max(0.0, r2)) / std::max(lam, 1e-12);
        }

        // Rounds of inexact inverse iteration + Rayleigh-Ritz; returns {rounds, residual}.
        std::pair<int, f64> resmooth(const CSRGraph &g, int budget, int cg_depth) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "resmooth");
            const std::vector<f64> deg = degrees(g);
            std::vector<f64> scratch(g.n);
            orthonormalise(0xA110u);
            int rounds = 0;
            f64 res = 0.0;
            for (int r = 0; r < budget; ++r) {
                for (auto &xk: X) {
                    scratch = xk;
                    cg_solve(g, deg, scratch, xk, cg_depth);
                    subtract_mean(xk);
                }
                orthonormalise(static_cast<u64>(r));
                rayleigh_ritz(g, deg);
                orthonormalise(static_cast<u64>(r) + 0x8000u);
                res = leading_residual(g, deg, scratch);
                ++rounds;
                if (res <= kResidualTol) break;
            }
            return {rounds, res};
        }

        // Exact inverse subspace iteration with a dense Cholesky factor of L + delta*I.
        void dense_solve(const CSRGraph &g) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "SpectralTracker", "dense_solve");
            const size_t n = g.n;
            const std::vector<f64> deg = degrees(g);
            const f64 maxdeg = deg.empty() ? 0.0 : *std::max_element(deg.begin(), deg.end());
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

            auto chol_solve = [&](std::vector<f64> &x) {
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
            };

            for (int r = 0; r < kCoarseRounds; ++r) {
                for (auto &xk: X) {
                    chol_solve(xk);
                    subtract_mean(xk);
                }
                orthonormalise(static_cast<u64>(r));
            }
            rayleigh_ritz(g, deg);
            orthonormalise(0x5EEDu);
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_SPECTRAL_TRACKER_H
