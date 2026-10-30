#ifndef HEIPROMAP_WAVER_CROWD_REFINEMENT_H
#define HEIPROMAP_WAVER_CROWD_REFINEMENT_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "waver_merge_refinement.h"

namespace HeiProMap {
    // Hyperparameters of the crowd refinement.
    class WaverCrowdRefinementConfiguration {
    public:
        explicit WaverCrowdRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_rounds = 5;      // rounds; stops earlier once a round commits no move
        f64 top_fraction = 0.2;  // the share of messages, most crowded first, that is worked on in a round
        f64 min_gain = 0.5;      // a move must lower the overlap sum by more than this
        bool verbose = false;
    };

    // Crowd refinement: shortens and separates the bounding boxes of the messages.
    //
    // The scheduler places a message against everything its bounding box can touch, so a message competes with every
    // other message whose box overlaps its own. A message is a pair of blocks that exchange halo values, and its box is
    // the rectangle spanned by the two blocks on the grid. The cover of a grid cell is the number of boxes that contain
    // it; the sum of the squared covers is the number of ordered pairs of boxes, weighted by the area they share, so
    // lowering it lowers both the cover (long boxes) and the crowd (boxes in busy regions).
    //
    // For the most crowded messages (mean cover inside their box), the sender group of either endpoint is moved to the
    // other endpoint or to one of its neighbors, if that lowers the sum of squared covers. The boxes that appear and
    // disappear by the move are found from the quotient rows, and the cover grid is updated as they are priced. A move
    // must not create more messages than it removes, and no block may get more neighbors than the busiest block
    // had at the start of the round, so the engine load does not pay for the shorter boxes.
    class WaverCrowdRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using Base = WaverMergeRefinement;

        template<typename DistanceOracleT>
        class CrowdEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::quotient;
            using E::collect_group; using E::measure; using E::apply; using E::targets_of;

            struct PairChange {
                partition_t p, q;   // p < q
                weight_t delta;
            };

        public:
            CrowdEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                        const WaverCrowdRefinementConfiguration &t_crowd,
                        const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, unused_config(), t_alpha_end, t_beta_end, t_alpha_link, t_beta_link),
                  crowd(t_crowd) {
                W = static_cast<int64_t>(t_oracle.current_W());
                H = (static_cast<int64_t>(this->k) + W - 1) / W;
                cover.assign(static_cast<size_t>(W * H), 0);
            }

            void run() {
                const partition_t k = this->k;
                for (partition_t a = 0; a < k; ++a) {
                    for (const auto &he: quotient.neighbors(a)) {
                        if (a < he.target) { change_box(a, he.target, +1); }
                    }
                }
                for (u64 round = 0; round < crowd.max_rounds; ++round) {
                    // the messages, most crowded first
                    std::vector<std::pair<f64, std::pair<partition_t, partition_t> > > messages;
                    size_t max_degree = 0;
                    f64 cover_sum = 0.0, crowd_sum = 0.0;
                    for (partition_t a = 0; a < k; ++a) {
                        max_degree = std::max(max_degree, quotient.neighbors(a).size());
                        for (const auto &he: quotient.neighbors(a)) {
                            if (a >= he.target) { continue; }
                            const f64 c = crowd_of(a, he.target);
                            messages.push_back({c, {a, he.target}});
                            crowd_sum += c;
                            cover_sum += static_cast<f64>(area_of(a, he.target));
                        }
                    }
                    if (messages.empty()) { return; }
                    std::sort(messages.begin(), messages.end(),
                              [](const auto &x, const auto &y) { return x.first > y.first; });
                    messages.resize(std::max<size_t>(1, static_cast<size_t>(crowd.top_fraction *
                                                                            static_cast<f64>(messages.size()))));
                    degree_cap = max_degree;

                    u64 tried = 0, kept = 0;
                    f64 gain = 0.0;
                    for (const auto &[c, ab]: messages) {
                        const auto [a, b] = ab;
                        if (quotient.weight(a, b) == 0) { continue; }
                        for (int side = 0; side < 2; ++side) {
                            f64 g_gain = 0.0;
                            const int r = try_move(side == 0 ? a : b, side == 0 ? b : a, g_gain);
                            if (r == 0) { continue; }
                            ++tried;
                            if (r > 0) {
                                ++kept;
                                gain += g_gain;
                                break;
                            }
                        }
                    }
                    if (crowd.verbose) {
                        std::cout << "[crowd] round " << round << ": " << messages.size() << " of "
                                  << static_cast<size_t>(messages.size() / std::max(crowd.top_fraction, 1e-9))
                                  << " messages, mean crowd " << crowd_sum / static_cast<f64>(messages.size() /
                                      std::max(crowd.top_fraction, 1e-9))
                                  << ", total area " << cover_sum << ", " << tried << " moves tried, " << kept
                                  << " kept, overlap gain " << gain << std::endl;
                    }
                    if (kept == 0) { break; }
                }
            }

        private:
            const WaverCrowdRefinementConfiguration &crowd;
            int64_t W = 1, H = 1;
            std::vector<int64_t> cover;
            size_t degree_cap = 0;

            static const WaverMergeRefinementConfiguration &unused_config() {
                static const WaverMergeRefinementConfiguration c("unused");
                return c;
            }

            struct Box {
                int64_t x0, x1, y0, y1;
            };

            Box box_of(const partition_t p, const partition_t q) const {
                const int64_t px = static_cast<int64_t>(p) % W, py = static_cast<int64_t>(p) / W;
                const int64_t qx = static_cast<int64_t>(q) % W, qy = static_cast<int64_t>(q) / W;
                return {std::min(px, qx), std::max(px, qx), std::min(py, qy), std::max(py, qy)};
            }

            int64_t area_of(const partition_t p, const partition_t q) const {
                const Box b = box_of(p, q);
                return (b.x1 - b.x0 + 1) * (b.y1 - b.y0 + 1);
            }

            // The mean cover inside the box of a message.
            f64 crowd_of(const partition_t p, const partition_t q) const {
                const Box b = box_of(p, q);
                int64_t sum = 0;
                for (int64_t y = b.y0; y <= b.y1; ++y) {
                    for (int64_t x = b.x0; x <= b.x1; ++x) { sum += cover[static_cast<size_t>(y * W + x)]; }
                }
                return static_cast<f64>(sum) / static_cast<f64>((b.x1 - b.x0 + 1) * (b.y1 - b.y0 + 1));
            }

            // Adds (+1) or removes (-1) the box of a message and returns the change of the sum of squared covers.
            int64_t change_box(const partition_t p, const partition_t q, const int sign) {
                const Box b = box_of(p, q);
                int64_t d = 0;
                for (int64_t y = b.y0; y <= b.y1; ++y) {
                    for (int64_t x = b.x0; x <= b.x1; ++x) {
                        int64_t &c = cover[static_cast<size_t>(y * W + x)];
                        if (sign > 0) {
                            d += 2 * c + 1;
                            ++c;
                        } else {
                            --c;
                            d -= 2 * c + 1;
                        }
                    }
                }
                return d;
            }

            // The change of the pair weights if the group moves from x to c.
            void pair_changes(const partition_t x, const partition_t c, const Base::Candidate &group,
                              std::vector<PairChange> &out) const {
                auto add = [&](partition_t p, partition_t q, weight_t d) {
                    if (p == q || d == 0) { return; }
                    if (p > q) { std::swap(p, q); }
                    out.push_back({p, q, d});
                };
                out.clear();
                for (const auto &[d, w]: group.ext) {
                    add(x, d, -w);
                    if (d != c) { add(c, d, w); }
                }
                add(x, c, group.weight_to_own);
                std::sort(out.begin(), out.end(), [](const PairChange &a, const PairChange &b) {
                    return a.p != b.p ? a.p < b.p : a.q < b.q;
                });
                size_t w = 0;
                for (size_t i = 0; i < out.size(); ++i) {
                    if (w > 0 && out[w - 1].p == out[i].p && out[w - 1].q == out[i].q) {
                        out[w - 1].delta += out[i].delta;
                    } else {
                        out[w++] = out[i];
                    }
                }
                out.resize(w);
            }

            // Prices the move of `group` from x to c: updates the cover grid for the boxes that appear and disappear,
            // and returns the change of the overlap sum. `ok` is false if the move would create more messages than it
            // removes or give some block too many neighbors. The grid is left as after the move; undo_boxes()
            // restores it.
            int64_t price(const partition_t x, const partition_t c, const Base::Candidate &group,
                      std::vector<std::pair<std::pair<partition_t, partition_t>, int> > &applied, bool &ok) {
                std::vector<PairChange> changes;
                pair_changes(x, c, group, changes);
                int64_t d_overlap = 0;
                int appear = 0, vanish = 0;
                ok = true;
                std::vector<std::pair<partition_t, int> > degree_change;
                for (const PairChange &ch: changes) {
                    const weight_t old_w = quotient.weight(ch.p, ch.q);
                    const weight_t new_w = old_w + ch.delta;
                    if (old_w == 0 && new_w > 0) {
                        ++appear;
                        d_overlap += change_box(ch.p, ch.q, +1);
                        applied.push_back({{ch.p, ch.q}, +1});
                        degree_change.push_back({ch.p, 1});
                        degree_change.push_back({ch.q, 1});
                    } else if (old_w > 0 && new_w <= 0) {
                        ++vanish;
                        d_overlap += change_box(ch.p, ch.q, -1);
                        applied.push_back({{ch.p, ch.q}, -1});
                        degree_change.push_back({ch.p, -1});
                        degree_change.push_back({ch.q, -1});
                    }
                }
                if (appear > vanish) { ok = false; }
                std::sort(degree_change.begin(), degree_change.end());
                for (size_t i = 0; i < degree_change.size() && ok;) {
                    size_t j = i;
                    int net = 0;
                    while (j < degree_change.size() && degree_change[j].first == degree_change[i].first) {
                        net += degree_change[j].second;
                        ++j;
                    }
                    if (net > 0 && quotient.neighbors(degree_change[i].first).size() + static_cast<size_t>(net) >
                                       degree_cap) {
                        ok = false;
                    }
                    i = j;
                }
                return d_overlap;
            }

            void undo_boxes(std::vector<std::pair<std::pair<partition_t, partition_t>, int> > &applied) {
                for (size_t i = applied.size(); i-- > 0;) {
                    change_box(applied[i].first.first, applied[i].first.second, -applied[i].second);
                }
                applied.clear();
            }

            // Moves the sender group of x towards y (to y or a neighbor of y) if that lowers the overlap sum.
            // Returns 0 if there is nothing to try, +1 if a move was kept, -1 if moves were tried and none kept.
            int try_move(const partition_t x, const partition_t y, f64 &gain) {
                Base::Candidate group;
                collect_group(x, y, group.group);
                if (group.group.empty()) { return 0; }
                measure(x, group);
                if (pm.bweights[x] - group.weight < pm.get_lmin(x)) { return 0; }
                std::vector<partition_t> targets;
                targets_of(x, y, targets);

                partition_t best = NO_ID;
                int64_t best_overlap = 0;
                bool any = false;
                std::vector<std::pair<std::pair<partition_t, partition_t>, int> > applied;
                for (const partition_t c: targets) {
                    if (pm.bweights[c] + group.weight > pm.get_lmax(c)) { continue; }
                    any = true;
                    bool ok = true;
                    const int64_t d = price(x, c, group, applied, ok);
                    undo_boxes(applied);
                    if (ok && static_cast<f64>(d) < -crowd.min_gain && d < best_overlap) {
                        best_overlap = d;
                        best = c;
                    }
                }
                if (!any) { return 0; }
                if (best == NO_ID) { return -1; }
                bool ok = true;
                price(x, best, group, applied, ok); // leaves the grid updated
                applied.clear();
                apply(group, x, best);
                gain = -static_cast<f64>(best_overlap);
                return 1;
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g, p_manager_t &p_manager, bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle, const WaverCrowdRefinementConfiguration &config,
                           const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverCrowdRefinement", "refine");
            CrowdEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                         alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_CROWD_REFINEMENT_H
