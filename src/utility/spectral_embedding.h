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

#ifndef HEIPROMAP_SPECTRAL_EMBEDDING_H
#define HEIPROMAP_SPECTRAL_EMBEDDING_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "profiler.h"

namespace HeiProMap {

    class SpectralEmbedding {
    public:
        // Compute the d smallest non-trivial eigenvectors of the normalized
        // graph Laplacian using inverse power iteration with deflation.
        // Returns eigenvectors as n x d matrix (row-major: result[i*d + dim]).
        static std::vector<f64> compute_laplacian_eigenvectors(
                const CSRGraph &g, size_t n_dims, u64 seed) {
            const size_t n = g.n;
            std::vector<f64> eigvecs(n * n_dims, 0.0);

            // Build dense adjacency-weighted Laplacian L = D - W
            // For small n (coarsest level), dense is fine
            std::vector<f64> L(n * n, 0.0);
            for (vertex_t u = 0; u < n; ++u) {
                f64 deg = 0.0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    const f64 w = static_cast<f64>(g.edge_w(j));
                    L[u * n + v] -= w;
                    deg += w;
                }
                L[u * n + u] = deg;
            }

            // Compute eigenvectors via repeated inverse iteration with shift
            // Use a small shift sigma to target smallest non-zero eigenvalues
            // (L - sigma*I) x = b, solve by Gaussian elimination

            // First, find all eigenvalues via Jacobi iteration on dense symmetric matrix
            // For n <= ~1000 this is practical
            std::vector<f64> A(L); // copy
            std::vector<f64> V(n * n, 0.0);
            for (size_t i = 0; i < n; ++i) V[i * n + i] = 1.0;

            jacobi_eigendecomposition(A, V, n);

            // Eigenvalues are on diagonal of A, eigenvectors in columns of V
            // Sort by eigenvalue ascending
            std::vector<size_t> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return A[a * n + a] < A[b * n + b];
            });

            // Skip the first eigenvector (constant, eigenvalue ~0), take next n_dims
            for (size_t d = 0; d < n_dims && d + 1 < n; ++d) {
                size_t col = order[d + 1];
                for (size_t i = 0; i < n; ++i) {
                    eigvecs[i * n_dims + d] = V[i * n + col];
                }
            }

            return eigvecs;
        }

        // Given a quotient graph and distance oracle, compute spectral embedding
        // of graph nodes and match to grid cells via min-cost assignment.
        // Returns permutation: perm[block] = cell.
        template<typename DistanceOracleT>
        static std::vector<partition_t> spectral_assignment(
                const CSRGraph &g,
                const PartitionManager &p_manager,
                const CSRGraph &quotient,
                const DistanceOracleT &d_oracle,
                u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("misc", "SpectralEmbedding", "spectral_assignment");
            const partition_t k = p_manager.get_k();
            const size_t W = static_cast<size_t>(d_oracle.current_W());

            // Compute 2D spectral embedding of quotient graph
            auto embedding = compute_laplacian_eigenvectors(quotient, 2, seed);

            // Grid cell coordinates (normalized to [0,1])
            const size_t H = (k + W - 1) / W;
            std::vector<f64> grid_x(k), grid_y(k);
            for (partition_t c = 0; c < k; ++c) {
                grid_x[c] = static_cast<f64>(c % W) / static_cast<f64>(std::max(W - 1, static_cast<size_t>(1)));
                grid_y[c] = static_cast<f64>(c / W) / static_cast<f64>(std::max(H - 1, static_cast<size_t>(1)));
            }

            // Normalize embedding to [0,1] per dimension
            for (size_t d = 0; d < 2; ++d) {
                f64 mn = embedding[d], mx = embedding[d];
                for (partition_t i = 0; i < k; ++i) {
                    mn = std::min(mn, embedding[i * 2 + d]);
                    mx = std::max(mx, embedding[i * 2 + d]);
                }
                f64 range = mx - mn;
                if (range < 1e-12) range = 1.0;
                for (partition_t i = 0; i < k; ++i) {
                    embedding[i * 2 + d] = (embedding[i * 2 + d] - mn) / range;
                }
            }

            // Try all 4 orientation variants (flip x, flip y) and both axis orderings
            // Pick the one with lowest matching cost
            std::vector<partition_t> best_perm(k);
            f64 best_cost = std::numeric_limits<f64>::max();

            for (int variant = 0; variant < 8; ++variant) {
                bool flip_x = (variant & 1) != 0;
                bool flip_y = (variant & 2) != 0;
                bool swap_axes = (variant & 4) != 0;

                std::vector<f64> ex(k), ey(k);
                for (partition_t i = 0; i < k; ++i) {
                    f64 vx = embedding[i * 2 + 0];
                    f64 vy = embedding[i * 2 + 1];
                    if (swap_axes) std::swap(vx, vy);
                    if (flip_x) vx = 1.0 - vx;
                    if (flip_y) vy = 1.0 - vy;
                    ex[i] = vx;
                    ey[i] = vy;
                }

                // Build cost matrix: cost[block][cell] = squared distance between
                // block's embedding and cell's grid position
                std::vector<std::vector<weight_t>> cost_matrix(k, std::vector<weight_t>(k));
                for (partition_t b = 0; b < k; ++b) {
                    for (partition_t c = 0; c < k; ++c) {
                        f64 dx = ex[b] - grid_x[c];
                        f64 dy = ey[b] - grid_y[c];
                        cost_matrix[b][c] = static_cast<weight_t>((dx * dx + dy * dy) * 1e6);
                    }
                }

                auto assignment = solve_hungarian(cost_matrix);

                // Compute actual QAP-like cost for this assignment
                f64 cost = 0.0;
                for (partition_t b = 0; b < k; ++b) {
                    for (partition_t c = 0; c < k; ++c) {
                        if (b == c) continue;
                        f64 dx = ex[b] - ex[c];
                        f64 dy = ey[b] - ey[c];
                        f64 embed_dist = std::sqrt(dx * dx + dy * dy);
                        f64 grid_dist = static_cast<f64>(d_oracle.get(assignment[b], assignment[c]));
                        cost += embed_dist * grid_dist;
                    }
                }

                if (cost < best_cost) {
                    best_cost = cost;
                    for (partition_t b = 0; b < k; ++b) {
                        best_perm[b] = static_cast<partition_t>(assignment[b]);
                    }
                }
            }

            return best_perm;
        }

    private:
        // Jacobi eigendecomposition for symmetric matrix.
        // A is overwritten with diagonal (eigenvalues), V gets eigenvectors as columns.
        static void jacobi_eigendecomposition(std::vector<f64> &A, std::vector<f64> &V, size_t n) {
            HEIPROMAP_PROFILE_SCOPE("misc", "SpectralEmbedding", "jacobi_eigendecomposition");
            const size_t max_sweeps = 100;
            const f64 eps = 1e-12;

            for (size_t sweep = 0; sweep < max_sweeps; ++sweep) {
                f64 off_diag = 0.0;
                for (size_t i = 0; i < n; ++i)
                    for (size_t j = i + 1; j < n; ++j)
                        off_diag += A[i * n + j] * A[i * n + j];

                if (off_diag < eps) break;

                for (size_t p = 0; p < n; ++p) {
                    for (size_t q = p + 1; q < n; ++q) {
                        if (std::abs(A[p * n + q]) < eps * 0.01) continue;

                        f64 tau = (A[q * n + q] - A[p * n + p]) / (2.0 * A[p * n + q]);
                        f64 t;
                        if (tau >= 0)
                            t = 1.0 / (tau + std::sqrt(1.0 + tau * tau));
                        else
                            t = -1.0 / (-tau + std::sqrt(1.0 + tau * tau));

                        f64 c = 1.0 / std::sqrt(1.0 + t * t);
                        f64 s = t * c;

                        f64 app = A[p * n + p];
                        f64 aqq = A[q * n + q];
                        f64 apq = A[p * n + q];

                        A[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                        A[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                        A[p * n + q] = 0.0;
                        A[q * n + p] = 0.0;

                        for (size_t r = 0; r < n; ++r) {
                            if (r == p || r == q) continue;
                            f64 arp = A[r * n + p];
                            f64 arq = A[r * n + q];
                            A[r * n + p] = c * arp - s * arq;
                            A[p * n + r] = A[r * n + p];
                            A[r * n + q] = s * arp + c * arq;
                            A[q * n + r] = A[r * n + q];
                        }

                        for (size_t r = 0; r < n; ++r) {
                            f64 vrp = V[r * n + p];
                            f64 vrq = V[r * n + q];
                            V[r * n + p] = c * vrp - s * vrq;
                            V[r * n + q] = s * vrp + c * vrq;
                        }
                    }
                }
            }
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_SPECTRAL_EMBEDDING_H
