#ifndef HEIPROMAP_SPECTRAL_PREPARTITIONER_H
#define HEIPROMAP_SPECTRAL_PREPARTITIONER_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "mapping.h"
#include "profiler.h"
#include "spectral_tracker.h"

namespace HeiProMap {

    // Partitions a graph into the k = W x H cells of a grid with exactly balanced, compact blocks, without an
    // external partitioner: a multilevel solve of the three smoothest Laplacian eigenvectors (heavy-edge
    // coarsening, dense solve at the coarsest level, prolongation and smoothing per level, all through
    // SpectralTracker), cut by recursive coordinate bisection with exact capacities. At every cut the
    // coordinate with the largest variance among the three is used, which keeps blocks compact for mesh-like
    // graphs (the three-coordinate spectral embedding of SpectralEmbeddingMapper).
    class SpectralPrepartitioner {
        static constexpr vertex_t kCoarsestN = 512;
        static constexpr f64 kStallRatio = 0.95;
        static constexpr vertex_t kNone = static_cast<vertex_t>(-1);

        // Greedy heavy-edge matching in vertex order (ties go to the lower vertex id); returns the fine -> coarse
        // vertex mapping. The deterministic order pairs grid-aligned neighbors, which keeps the coarse levels
        // regular. A random visiting order gave 18% more messages and 20% more halo volume on the 10M-row Poisson test.
        static Mapping heavy_edge_matching(const CSRGraph &g) {
            const vertex_t n = g.n;

            std::vector<vertex_t> match(n, kNone);
            for (vertex_t u = 0; u < n; ++u) {
                if (match[u] != kNone) { continue; }
                vertex_t best = kNone;
                weight_t best_w = -1;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    if (v == u || match[v] != kNone) { continue; }
                    const weight_t w = g.edge_w(j);
                    if (w > best_w || (w == best_w && v < best)) {
                        best_w = w;
                        best = v;
                    }
                }
                if (best != kNone) {
                    match[u] = best;
                    match[best] = u;
                } else {
                    match[u] = u;
                }
            }

            Mapping mapping;
            mapping.initialize(n);
            vertex_t coarse_n = 0;
            for (vertex_t u = 0; u < n; ++u) {
                if (match[u] >= u) {
                    mapping.set(u, coarse_n);
                    if (match[u] != u) { mapping.set(match[u], coarse_n); }
                    ++coarse_n;
                }
            }
            mapping.set_coarse_n(coarse_n);
            return mapping;
        }

        // Recursive coordinate bisection of the vertices onto the W x H grid with exact capacities (a side
        // receives members * cells_in_side / cells_in_rectangle vertices, rounded down). Block id = row * W + col.
        static std::vector<partition_t> bisect(const std::vector<SpectralTracker::Vec> &X, const vertex_t n,
                                               const size_t W, const size_t H, const size_t dims) {
            HEIPROMAP_PROFILE_SCOPE("partition", "SpectralPrepartitioner", "bisect");
            struct Node {
                size_t lo, hi, col0, row0, w, h;
            };
            std::vector<partition_t> out(n, 0);
            std::vector<vertex_t> perm(n);
            std::iota(perm.begin(), perm.end(), static_cast<vertex_t>(0));

            std::vector<Node> stack{{0, static_cast<size_t>(n), 0, 0, W, H}};
            while (!stack.empty()) {
                const Node nd = stack.back();
                stack.pop_back();

                if (nd.w * nd.h == 1) {
                    const auto cell = static_cast<partition_t>(nd.row0 * W + nd.col0);
                    for (size_t i = nd.lo; i < nd.hi; ++i) { out[perm[i]] = cell; }
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

                // Cut along the coordinate with the largest spread among the members.
                const size_t m = nd.hi - nd.lo;
                size_t coord = 0;
                f64 best_var = -1.0;
                // dims 0: column cuts use the smoothest vector, row cuts the second (the longer side of the grid gets
                // the smoothest), as in SpectralTracker::place, so that the blocks line up with the rows and columns.
                if (dims == 0) { coord = (split_cols == (W >= H)) ? 0 : 1; }
                for (size_t c = 0; dims > 0 && c < std::min(dims, X.size()); ++c) {
                    f64 sum = 0.0, sumsq = 0.0;
                    for (size_t i = nd.lo; i < nd.hi; ++i) {
                        const f64 x = X[c][perm[i]];
                        sum += x;
                        sumsq += x * x;
                    }
                    const f64 var = sumsq - sum * sum / static_cast<f64>(m);
                    if (var > best_var) {
                        best_var = var;
                        coord = c;
                    }
                }

                const size_t cells = nd.w * nd.h;
                const size_t cells_a = a.w * a.h;
                const size_t count_a = static_cast<size_t>(static_cast<u64>(m) * cells_a / cells);
                const SpectralTracker::Vec &key = X[coord];
                std::nth_element(perm.begin() + nd.lo, perm.begin() + nd.lo + count_a, perm.begin() + nd.hi,
                                 [&](vertex_t x, vertex_t y) {
                                     return key[x] != key[y] ? key[x] < key[y] : x < y;
                                 });

                a.lo = nd.lo;
                a.hi = nd.lo + count_a;
                b.lo = a.hi;
                b.hi = nd.hi;
                stack.push_back(b);
                stack.push_back(a);
            }
            return out;
        }

        static f64 cut_weight(const CSRGraph &g, const std::vector<partition_t> &part) {
            f64 cut = 0.0;
            for (vertex_t u = 0; u < g.n; ++u) {
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    if (part[u] != part[g.edges_v[j]]) { cut += static_cast<f64>(g.edge_w(j)); }
                }
            }
            return cut;
        }

    public:
        static constexpr size_t kAuto = static_cast<size_t>(-1);

        // Returns part[u] in [0, W * H) for every vertex of g. Requires g.n >= W * H.
        static std::vector<partition_t> partition(const CSRGraph &g, const size_t W, const size_t H,
                                                  const u64 seed, const u64 threads,
                                                  const SpectralSmoothing &smoothing = SpectralSmoothing(),
                                                  const u64 rotations = 1, const size_t dims = 3) {
            HEIPROMAP_PROFILE_SCOPE("partition", "SpectralPrepartitioner", "partition");

            // Heavy-edge coarsening hierarchy. Reserved so the graph references stay valid.
            std::vector<CSRGraph> coarse;
            std::vector<Mapping> maps; // maps[i]: level i (fine) -> level i + 1 (coarse); level 0 is g
            coarse.reserve(64);
            maps.reserve(64);
            auto level = [&](const size_t i) -> const CSRGraph & { return i == 0 ? g : coarse[i - 1]; };
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "SpectralPrepartitioner", "coarsen");
                while (level(maps.size()).n > kCoarsestN && coarse.size() < 63) {
                    const CSRGraph &fine = level(maps.size());
                    Mapping mapping = heavy_edge_matching(fine);
                    if (static_cast<f64>(mapping.get_coarse_n()) > kStallRatio * static_cast<f64>(fine.n)) { break; }
                    maps.push_back(std::move(mapping));
                    coarse.emplace_back();
                    coarse.back().contract(fine, maps.back(), threads);
                }
            }

            SpectralTracker tracker;
            SpectralSmoothing tracker_smoothing = smoothing;
            tracker_smoothing.threads = static_cast<int>(threads);
            tracker.set_smoothing(tracker_smoothing);
            tracker.init(level(maps.size()), seed);
            for (size_t i = maps.size(); i-- > 0;) { tracker.prolong(maps[i], level(i)); }

            // dims kAuto: the cuts by largest spread of all three vectors (3D-like graphs) and the cuts along the rows
            // and columns by the two smoothest (2D-like graphs) compete, the smaller cut wins.
            std::vector<partition_t> best = bisect(tracker.basis(), g.n, W, H, dims == kAuto ? 3 : dims);
            if (dims == kAuto) {
                std::vector<partition_t> aligned = bisect(tracker.basis(), g.n, W, H, 0);
                if (cut_weight(g, aligned) < cut_weight(g, best)) { best = std::move(aligned); }
            }
            if (rotations <= 1) { return best; }

            // The two smoothest vectors of a (nearly) square domain are an arbitrary turn of the axis-aligned ones, and
            // cutting along a turned pair gives slanted blocks: try turns of that pair over a quarter turn and keep the
            // partition with the smallest cut.
            f64 best_cut = cut_weight(g, best);
            const std::vector<SpectralTracker::Vec> &X = tracker.basis();
            std::vector<SpectralTracker::Vec> turned(X.size());
            for (auto &t: turned) { t.initialize(g.n); }
            for (size_t c = 2; c < X.size(); ++c) { for (vertex_t u = 0; u < g.n; ++u) { turned[c][u] = X[c][u]; } }
            for (u64 step = 1; step < rotations; ++step) {
                const f64 theta = 0.5 * 3.14159265358979323846 * static_cast<f64>(step) / static_cast<f64>(rotations);
                const f64 cs = std::cos(theta), sn = std::sin(theta);
                for (vertex_t u = 0; u < g.n; ++u) {
                    turned[0][u] = cs * X[0][u] + sn * X[1][u];
                    turned[1][u] = -sn * X[0][u] + cs * X[1][u];
                }
                std::vector<partition_t> candidate = bisect(turned, g.n, W, H, dims == kAuto ? 3 : dims);
                const f64 c = cut_weight(g, candidate);
                if (c < best_cut) {
                    best_cut = c;
                    best = std::move(candidate);
                }
            }
            return best;
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_SPECTRAL_PREPARTITIONER_H
