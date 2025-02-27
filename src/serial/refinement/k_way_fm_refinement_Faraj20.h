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
        u64 max_iteration = 1; // how many iterations to run the algorithm at most
        f64 alpha         = 1000.0;
        f64 beta          = 1.0;
    };

    class KWayFMMove {
    public:
        vertex_t    u;
        partition_t u_id;
        partition_t to_move_id;
        s64         qap_delta;

        KWayFMMove(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            u          = t_u;
            u_id       = t_u_id;
            to_move_id = t_to_move;
            qap_delta  = t_qap_delta;
        }

        bool operator>(const KWayFMMove &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const KWayFMMove &m) const {
            return qap_delta < m.qap_delta;
        }
    };

    class KWayFMMoves {
    private:
        vertex_t                m_u;
        s64                     max_qap_delta = std::numeric_limits<s64>::min();
        std::vector<KWayFMMove> moves;

    public:
        KWayFMMoves(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            m_u           = t_u;
            max_qap_delta = t_qap_delta;
            push(t_u, t_u_id, t_to_move, t_qap_delta);
        }

        void push(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            max_qap_delta = std::max(max_qap_delta, t_qap_delta);
            for (auto &move: moves) {
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

        void swap(KWayFMMoves &other) noexcept {
            std::swap(m_u, other.m_u);
            std::swap(max_qap_delta, other.max_qap_delta);
            moves.swap(other.moves);
        }

        bool operator<(const KWayFMMoves &m) const {
            return get_max_qap_delta() < m.get_max_qap_delta();
        }

        bool operator>(const KWayFMMoves &m) const {
            return get_max_qap_delta() > m.get_max_qap_delta();
        }

        bool operator==(const KWayFMMoves &m) const {
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
        vertex_t                 m_n = 0;
        std::vector<KWayFMMoves> m_heap; // The heap array
        std::vector<size_t>      m_indices; // Mapping of keys to heap indices

        u64              m_iteration = 1;
        std::vector<u64> m_iteration_counter;

    public:
        size_t push_operations = 0;

        KWayFMPriorityQueue() = default;

        explicit KWayFMPriorityQueue(const vertex_t t_n) {
            m_n = t_n;

            m_heap.reserve(m_n);
            m_indices.resize(m_n);
            m_iteration_counter.resize(m_n);
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

            m_heap[m_indices[t_u]].sort();
            bubble_up(m_indices[t_u]);
            bubble_down(m_indices[t_u]);
        }

        KWayFMMove top() const {
            ASSERT(!m_heap.empty());
            ASSERT(!m_heap[0].empty());
            return m_heap[0].top();
        }

        void pop() {
            ASSERT(!empty());

            swap(0, m_heap.size() - 1);
            m_indices[m_heap.back().get_u()] = HEAP_TOMBSTONE;
            m_heap.pop_back();
            if (!m_heap.empty()) {
                bubble_down(0);
            }
        }

        bool empty() const {
            return m_heap.empty();
        }

        void clear() {
            push_operations = 0;
            m_heap.clear();
            m_iteration += 1;
        }

    private:
        bool entry_exists(const vertex_t t_u) const {
            ASSERT(t_u < m_n);
            return m_iteration_counter[t_u] == m_iteration && m_indices[t_u] != HEAP_TOMBSTONE;
        }

        void push_new(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            ASSERT(!entry_exists(t_u));
            m_heap.emplace_back(t_u, t_u_id, t_to_move, t_qap_delta);
            m_indices[t_u]           = m_heap.size() - 1;
            m_iteration_counter[t_u] = m_iteration;
            bubble_up(m_heap.size() - 1);
        }

        void update(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const s64 t_qap_delta) {
            ASSERT(entry_exists(t_u));
            m_heap[m_indices[t_u]].push(t_u, t_u_id, t_to_move, t_qap_delta);
            bubble_up(m_indices[t_u]);
            bubble_down(m_indices[t_u]);
        }

        // Bubbles up the element at the given index to restore the heap property
        void bubble_up(size_t index) {
            while (index > 0) {
                size_t parent_index = (index - 1) / 2;
                if (m_heap[index] <= m_heap[parent_index]) break;

                swap(index, parent_index);
                index = parent_index;
            }
        }

        // Bubbles down the element at the given index to restore the heap property
        void bubble_down(size_t index) {
            size_t size = m_heap.size();
            ASSERT(index < size);

            while (true) {
                size_t left_child_index  = 2 * index + 1;
                size_t right_child_index = 2 * index + 2;
                size_t largest_index     = index;

                if (left_child_index < size && m_heap[left_child_index] > m_heap[largest_index]) {
                    largest_index = left_child_index;
                }
                if (right_child_index < size && m_heap[right_child_index] > m_heap[largest_index]) {
                    largest_index = right_child_index;
                }
                if (largest_index == index) break;

                swap(index, largest_index);
                index = largest_index;
            }
        }

        void swap(const size_t i, const size_t j) {
            m_heap[i].swap(m_heap[j]);
            m_indices[m_heap[i].get_u()] = i;
            m_indices[m_heap[j].get_u()] = j;
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
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        u64                      m_seed = 0;

        std::vector<u32> vertex_used;
        u32              vertex_mark = 0;

        std::vector<u32> block_used;
        u32              block_marker = 0;

        // indexed max heaps
        KWayFMPriorityQueue queue;

        std::priority_queue<KWayFMMove> prio_queue;

        std::mt19937                          gen;
        std::uniform_real_distribution<float> dis;

    public:
        KWayFMRefinementFaraj20() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_seed) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_seed      = t_seed;

            vertex_used.resize(t_n, 0);
            block_used.resize(t_n, 0);

            queue = KWayFMPriorityQueue(t_n);

            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] KWayFMRefinementFaraj20Configuration &config,
                    [[maybe_unused]] TSerialGraph &g,
                    [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager &bv_manager,
                    [[maybe_unused]] TSerialPartitionManager &p_manager,
                    [[maybe_unused]] TSerialDistanceOracle &d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph &q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            std::vector<KWayFMMove> moves;

            config.beta = std::log(av_manager.get_n_active());

            // std::cout << "alpha = " << config.alpha << std::endl;
            // std::cout << "beta  = " << config.beta << std::endl;

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                auto sp    = std::chrono::high_resolution_clock::now();

                queue.clear();
                vertex_mark += 1;
                prio_queue = std::priority_queue<KWayFMMove>();

                // insert all boundary vertices
                for (vertex_t u: bv_manager) {
                    partition_t u_id     = p_manager[u];
                    weight_t    u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    block_marker += 1;
                    bool one_id_is_valid = false;
                    for (const auto [v, w]: g[u]) {
                        partition_t v_id = p_manager[v];
                        if (v_id == u_id) { continue; }
                        if (block_used[v_id] == block_marker) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        // queue.push(u, u_id, v_id, qap_delta);
                        prio_queue.emplace(u, u_id, v_id, qap_delta);
                        one_id_is_valid = true;

                        block_used[v_id] = block_marker;
                    }
                    if (one_id_is_valid) {
                        // queue.sort(u);
                    }
                }

                moves.clear();
                size_t best_idx      = 0;
                s64    curr_qap_gain = 0;
                s64    max_qap_gain  = 0;

                // u64 steps_since_last_improvement = 0;
                // f64 qap_gain_mean                = 0.0;
                // f64 qap_gain_var                 = 1.0;

                while (!prio_queue.empty()) {
                    const KWayFMMove move = prio_queue.top();
                    prio_queue.pop();

                    vertex_t vertex = move.u;
                    if (vertex_used[vertex] == vertex_mark) { continue; }

                    partition_t move_id = move.to_move_id;
                    if (!is_connected_to(g, p_manager, vertex, move_id)) { continue; }

                    partition_t vertex_id     = p_manager[vertex];
                    weight_t    vertex_weight = g.get_weight(vertex);
                    if (p_manager.get_bweight(move_id) + vertex_weight > m_lmax) { continue; }

                    s64 curr_qap_delta = get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle);
                    if (curr_qap_delta != move.qap_delta) { continue; }

                    moves.push_back(move);
                    curr_qap_gain += move.qap_delta;
                    if (curr_qap_gain >= max_qap_gain) {
                        best_idx     = moves.size();
                        max_qap_gain = curr_qap_gain;

                        // steps_since_last_improvement = 0;
                        // qap_gain_mean                = 0.0;
                        // qap_gain_var                 = 1.0;
                    }

                    // make move in structures
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    vertex_used[vertex] = vertex_mark;

                    // steps_since_last_improvement += 1;
                    // f64 new_qap_gain_mean = qap_gain_mean + ((f64) move.qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                    // f64 new_qap_gain_var  = (qap_gain_var + ((f64) move.qap_delta - qap_gain_mean) * ((f64) move.qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                    // qap_gain_mean = new_qap_gain_mean;
                    // qap_gain_var  = new_qap_gain_var;

                    /*
                    if (steps_since_last_improvement > 3 && (f64)steps_since_last_improvement * qap_gain_mean * qap_gain_mean > config.alpha * qap_gain_var + config.beta) {
                        std::cout << "Stop on random walk: " << steps_since_last_improvement << " " << qap_gain_mean << " " << qap_gain_var << std::endl;
                        break;
                    }
                    */


                    // we have to push or update the neighbors that were not moved already
                    for (const auto [neighbor, _]: g[vertex]) {
                        if (vertex_used[neighbor] == vertex_mark) { continue; }
                        if (!is_boundary(g, p_manager, neighbor)) { continue; }

                        partition_t neighbor_id     = p_manager[neighbor];
                        weight_t    neighbor_weight = g.get_weight(neighbor);

                        block_marker += 1;
                        bool one_id_is_valid = false;
                        for (const auto [v, w]: g[neighbor]) {
                            partition_t v_id = p_manager[v];
                            if (v_id == neighbor_id) { continue; }
                            if (block_used[v_id] == block_marker) { continue; }
                            if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                            // queue.push(neighbor, neighbor_id, v_id, qap_delta);
                            prio_queue.emplace(neighbor, neighbor_id, v_id, qap_delta);
                            one_id_is_valid = true;
                            block_used[v_id] = block_marker;
                        }
                        if (one_id_is_valid) {
                            queue.sort(neighbor);
                        }
                    }
                }

                // revert all moves in partitioning manager
                for (size_t i = 0; i < moves.size(); i++) {
                    vertex_t    vertex        = moves[moves.size() - 1 - i].u;
                    weight_t    vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id     = moves[moves.size() - 1 - i].to_move_id;
                    partition_t move_id       = moves[moves.size() - 1 - i].u_id;

                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }

                // make all moves to best index
                for (size_t i = 0; i < best_idx; ++i) {
                    vertex_t    vertex        = moves[i].u;
                    weight_t    vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id     = moves[i].u_id;
                    partition_t move_id       = moves[i].to_move_id;

                    bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                    q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }

                auto ep = std::chrono::high_resolution_clock::now();
                std::cout << "iteration: " << iteration << " best_idx: " << best_idx << " new_gain: " << max_qap_gain << " push ops: " << queue.push_operations << " max moves: " << moves.size() << " time: " << get_seconds(sp, ep) << std::endl;
            }
        }
    };
}


#endif //HEIPROMAP_K_WAY_FM_REFINEMENT_FARAJ20_H
