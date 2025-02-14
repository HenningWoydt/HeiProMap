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

#ifndef HEIPROMAP_K_WAY_FM_REFINEMENT_FARAJ20_H
#define HEIPROMAP_K_WAY_FM_REFINEMENT_FARAJ20_H

#include <algorithm>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct KWayFMRefinementFaraj20Configuration {
        u64 max_iteration = 2; // how many iterations to run the algorithm at most
    };

    class KWayFMMove {
    public:
        vertex_t u;
        partition_t u_id;
        partition_t to_move_id;
        s64 qap_delta;

        KWayFMMove(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            u          = t_u;
            u_id       = t_u_id;
            to_move_id = t_to_move;
            qap_delta  = t_qap_delta;
        }

        bool operator>(const KWayFMMove& m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const KWayFMMove& m) const {
            return qap_delta < m.qap_delta;
        }
    };

    class KWayFMMoves {
    private:
        vertex_t m_u;
        s64 max_qap_delta = std::numeric_limits<s64>::min();
        std::vector<KWayFMMove> moves;

    public:
        KWayFMMoves(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            m_u           = t_u;
            max_qap_delta = t_qap_delta;
            push(t_u, t_u_id, t_to_move, t_qap_delta);
        }

        void push(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            max_qap_delta = std::max(max_qap_delta, t_qap_delta);
            for (auto& move : moves) {
                if (move.to_move_id == t_to_move) {
                    move.qap_delta = t_qap_delta;
                    return;
                }
            }
            moves.emplace_back(t_u, t_u_id, t_to_move, t_qap_delta);
        }

        void sort() {
            std::sort(moves.begin(), moves.end());
        }

        KWayFMMove top() const {
            return moves.back();
        }

        bool empty() const {
            return moves.empty();
        }

        vertex_t get_u() const {
            return m_u;
        }

        s64 get_max_qap_delta() const {
            return max_qap_delta;
        }

        bool operator<(const KWayFMMoves& m) const {
            return get_max_qap_delta() < m.get_max_qap_delta();
        }

        bool operator>(const KWayFMMoves& m) const {
            return get_max_qap_delta() > m.get_max_qap_delta();
        }

        bool operator==(const KWayFMMoves& m) const {
            return get_max_qap_delta() == m.get_max_qap_delta();
        }

        bool operator<=(std::vector<KWayFMMoves>::const_reference value) const {
            return get_max_qap_delta() <= value.get_max_qap_delta();
        }

        bool operator>=(std::vector<KWayFMMoves>::const_reference value) const {
            return get_max_qap_delta() >= value.get_max_qap_delta();
        }
    };

    class KWayFMPriorityQueue {
        vertex_t m_n = 0;
        std::vector<KWayFMMoves> heap; // The heap array
        std::vector<size_t> indices; // Mapping of keys to heap indices

    public:
        size_t push_operations = 0;
        KWayFMPriorityQueue() = default;

        explicit KWayFMPriorityQueue(const vertex_t t_n) {
            m_n = t_n;

            heap.reserve(m_n);
            indices.resize(m_n, HEAP_TOMBSTONE);
        }

        void push(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            push_operations++;
            if (entry_exists(t_u)) {
                update(t_u, t_u_id, t_to_move, t_qap_delta);
            } else {
                push_new(t_u, t_u_id, t_to_move, t_qap_delta);
            }
        }

        void sort(const vertex_t t_u) {
            ASSERT(entry_exists(t_u));

            heap[indices[t_u]].sort();
            bubble_up(indices[t_u]);
            bubble_down(indices[t_u]);
        }

        KWayFMMove top() const {
            ASSERT(!heap.empty());
            ASSERT(!heap[0].empty());
            return heap[0].top();
        }

        void pop() {
            ASSERT(!empty());

            swap(0, heap.size() - 1);
            indices[heap.back().get_u()] = HEAP_TOMBSTONE;
            heap.pop_back();
            if (!heap.empty()) {
                bubble_down(0);
            }
        }

        bool empty() const {
            return heap.empty();
        }

        void clear() {
            push_operations = 0;
            heap.clear();
            std::fill(indices.begin(), indices.end(), HEAP_TOMBSTONE);
        }

    private:
        bool entry_exists(const vertex_t t_u) const {
            ASSERT(t_u < m_n);
            return indices[t_u] != HEAP_TOMBSTONE;
        }

        void push_new(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            ASSERT(!entry_exists(t_u));
            heap.emplace_back(t_u, t_u_id, t_to_move, t_qap_delta);
            indices[t_u] = heap.size() - 1;
            bubble_up(heap.size() - 1);
        }

        void update(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            ASSERT(entry_exists(t_u));
            heap[indices[t_u]].push(t_u, t_u_id, t_to_move, t_qap_delta);
            bubble_up(indices[t_u]);
            bubble_down(indices[t_u]);
        }

        // Bubbles up the element at the given index to restore the heap property
        void bubble_up(size_t index) {
            while (index > 0) {
                size_t parent_index = (index - 1) / 2;
                if (heap[index] <= heap[parent_index]) break;

                swap(index, parent_index);
                index = parent_index;
            }
        }

        // Bubbles down the element at the given index to restore the heap property
        void bubble_down(size_t index) {
            size_t size = heap.size();
            ASSERT(index < size);

            while (true) {
                size_t left_child_index  = 2 * index + 1;
                size_t right_child_index = 2 * index + 2;
                size_t largest_index     = index;

                if (left_child_index < size && heap[left_child_index] > heap[largest_index]) {
                    largest_index = left_child_index;
                }
                if (right_child_index < size && heap[right_child_index] > heap[largest_index]) {
                    largest_index = right_child_index;
                }
                if (largest_index == index) break;

                swap(index, largest_index);
                index = largest_index;
            }
        }

        void swap(const size_t i, const size_t j) {
            std::swap(heap[i], heap[j]);
            indices[heap[i].get_u()] = i;
            indices[heap[j].get_u()] = j;
        }
    };

    /**
     * Executes K-Way FM Refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania, Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     * > Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2020.
     *
     * This includes stopping based on a random walk described in
     * > Sanders, P., Schulz, C. (2011).
     * > Engineering Multilevel Graph Partitioning Algorithms.
     * > In: Demetrescu, C., Halldórsson, M.M. (eds) Algorithms – ESA 2011. ESA 2011. Lecture Notes in Computer Science, vol 6942. Springer, Berlin, Heidelberg
     *
     */
    class KWayFMRefinementFaraj20 final : public ISerialRefiner {
    private:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        std::vector<s32> used;
        s32 mark = -1;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

        // indexed max heaps
        KWayFMPriorityQueue queue;
        // KWayFMPriorityQueueNew queue_new;

    public:
        KWayFMRefinementFaraj20() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t n,
                        std::vector<partition_t>& t_hierarchy,
                        std::vector<weight_t>& t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);

            queue     = KWayFMPriorityQueue(n);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] KWayFMRefinementFaraj20Configuration& config,
                    [[maybe_unused]] TSerialGraph& g,
                    [[maybe_unused]] TSerialActiveVertexManager& av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager& bv_manager,
                    [[maybe_unused]] TSerialPartitionManager& p_manager,
                    [[maybe_unused]] TSerialDistanceOracle& d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph& q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            u32 counter = 0;
            std::vector<u32> found_ids_mark(k, counter);

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                auto sp = std::chrono::high_resolution_clock::now();
                queue.clear();
                mark += 1;

                // insert all boundary vertices
                for (vertex_t u : bv_manager) {
                    partition_t u_id  = p_manager[u];
                    weight_t u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    counter += 1;
                    bool one_id_is_valid = false;
                    for (const auto& [v, w] : g[u]) {
                        // for (size_t idx = 0; idx < g.size(u); ++idx) {
                        // const vertex_t v = g.neighbor(u, idx);
                        // const weight_t w = g.get_weight(u, idx);
                        partition_t v_id = p_manager[v];
                        if (v_id != u_id && found_ids_mark[v_id] != counter && p_manager.get_bweight(v_id) + u_weight <= lmax) {
                            found_ids_mark[v_id] = counter;
                            s64 qap_delta        = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                            if (qap_delta <= 0) {
                                // at first skip negative gains
                                continue;
                            }

                            queue.push(u, u_id, v_id, qap_delta);
                            one_id_is_valid = true;
                        }
                    }
                    if (one_id_is_valid) {
                        queue.sort(u);
                    }
                }

                std::vector<KWayFMMove> moves;
                size_t best_idx   = 0;
                s64 curr_qap_gain = 0;
                s64 max_qap_gain  = 0;

                while (!queue.empty()) {
                    const KWayFMMove move     = queue.top();
                    queue.pop();

                    vertex_t vertex        = move.u;
                    weight_t vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id  = p_manager[vertex];
                    partition_t move_id    = move.to_move_id;

                    if (vertex_id != move.u_id || used[vertex] == mark || p_manager.get_bweight(move_id) + vertex_weight > lmax) {
                        // if vertex_id and old vertex_id don't match
                        // vertex was already used
                        // moving would overload
                        continue;
                    }

                    used[vertex] = mark;

                    // make move in structures
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                    moves.push_back(move);
                    curr_qap_gain += move.qap_delta;
                    if (curr_qap_gain > max_qap_gain && !p_manager.is_overloaded()) {
                        best_idx     = moves.size();
                        max_qap_gain = curr_qap_gain;
                    }

                    // we have to push or update the neighbors that were not moved already
                    for (const auto& [neighbor, w] : g[vertex]) {
                        // for (size_t idx = 0; idx < g.size(vertex); ++idx) {
                        // const vertex_t neighbor = g.neighbor(vertex, idx);
                        // const weight_t w = g.get_weight(vertex, idx);
                        if (used[neighbor] == mark || !is_boundary(g, p_manager, vertex)) {
                            continue;
                        }
                        partition_t neighbor_id = p_manager[neighbor];

                        counter += 1;
                        bool one_id_is_valid = false;
                        for (const auto& [v, w] : g[neighbor]) {
                            partition_t v_id = p_manager[v];
                            if (v_id != neighbor_id && found_ids_mark[v_id] != counter) {
                                found_ids_mark[v_id] = counter;
                                s64 qap_delta        = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);

                                queue.push(neighbor, neighbor_id, v_id, qap_delta);
                                one_id_is_valid = true;
                            }
                        }
                        if (one_id_is_valid) {
                            queue.sort(neighbor);
                        }
                    }
                }

                // revert all moves in partitioning manager
                for (size_t i = 0; i < moves.size(); i++) {
                    vertex_t vertex        = moves[moves.size() - 1 - i].u;
                    weight_t vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id  = moves[moves.size() - 1 - i].to_move_id;
                    partition_t move_id    = moves[moves.size() - 1 - i].u_id;

                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }

                // make all moves to best index
                for (size_t i = 0; i < best_idx; ++i) {
                    vertex_t vertex        = moves[i].u;
                    weight_t vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id  = moves[i].u_id;
                    partition_t move_id    = moves[i].to_move_id;

                    bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                    q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }

                auto ep = std::chrono::high_resolution_clock::now();
                std::cout << "iteration: " << iteration << " " << best_idx << " " << max_qap_gain << " push ops: " << queue.push_operations << " max moves: " << moves.size() << " time: " << get_seconds(sp, ep) << std::endl;
            }
        }
    };
}


#endif //HEIPROMAP_K_WAY_FM_REFINEMENT_FARAJ20_H
