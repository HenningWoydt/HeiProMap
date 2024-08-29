#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_SOLVER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_SOLVER_H

#include "../../serial/datastructures/graph.h"
#include "../../serial/datastructures/active_vertex_manager.h"
#include "../../serial/datastructures/partition_manager.h"
#include "../../serial/datastructures/boundary_vertex_manager.h"
#include "../../serial/coarsening/heavy_edge_matcher.h"
#include "../../serial/partitioning/kaffpa_partitioner.h"
#include "../../serial/refinement/label_propagation_refinement.h"
#include "../../serial/refinement/identity_refinement.h"
#include "../../serial/datastructures/statistic_collector.h"
#include "parallel_boundary_vertex_manager.h"
#include "parallel_partition_manager.h"
#include "parallel_active_vertex_manager.h"
#include "parallel_graph.h"
#include "../coarsening/parallel_heavy_edge_matcher.h"
#include "../partitioning/parallel_kaffpa_partitioner.h"
#include "../refinement/parallel_label_propagation_refinement.h"
#include "parallel_distance_oracle.h"
#include "parallel_quotient_graph.h"

namespace HeiProMap {
    class ParallelSolver {
    private:
        // main structures
        ParallelGraph m_g;
        ParallelActiveVertexManager m_av_manager;
        ParallelPartitionManager m_p_manager;
        ParallelBoundaryVertexManager m_bv_manager;
        ParallelQuotientGraph m_qgraph;

        // distance
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k;
        ParallelDistanceOracle m_d_oracle;

        // balance
        weight_t m_lmax = 0;

        // threads
        u64 m_n_threads = 1;

        // multilevel
        vertex_t m_threshold;

        // matching
        // GreedyEdgeMatcher gem;
        ParallelHeavyEdgeMatcher m_he_matcher;
        // SimpleClustering sc;

        std::vector<std::vector<EdgeUV>> m_matches;

        // partitioning
        ParallelKaffpaPartitioner m_kaffpa_partitioner;

        // refinement
        ParallelLabelPropagationRefinement m_lp_refine;
        // IdentityRefinement m_i_refine;
        // QuotientGraphRefinement qgr;

        // statistics
        StatisticCollector m_stat_collect;

    public:
        ParallelSolver(std::string &t_graph_in,
                       std::vector<partition_t> &t_hierarchy,
                       std::vector<weight_t> &t_distance,
                       f64 t_imbalance,
                       u64 t_n_threads) {
            auto sp_graph_io = std::chrono::high_resolution_clock::now();
            m_g.initialize(t_graph_in, t_n_threads);
            auto ep_graph_io = std::chrono::high_resolution_clock::now();

            auto sp_io = std::chrono::high_resolution_clock::now();

            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_k = prod<partition_t>(m_hierarchy);

            // manager
            m_av_manager.initialize(&m_g, t_n_threads);
            m_p_manager.initialize(&m_g, &m_av_manager, m_k, t_n_threads);
            m_bv_manager.initialize(&m_g, &m_av_manager, &m_p_manager, m_k, t_n_threads);
            m_qgraph.initialize(m_k, t_n_threads);
            HEAVYASSERT(assert_state_pre_partitioning(m_g, m_av_manager));

            // distance
            m_d_oracle.initialize(m_hierarchy, m_distance, t_n_threads);

            // balance
            m_lmax = ceil((1.0 + t_imbalance) * ((f64) m_g.get_weight() / (f64) m_k));

            // threads
            m_n_threads = t_n_threads;

            // multilevel
            m_threshold = m_g.get_n() / 500;

            // matching
            // gem.initialize(&g);
            m_he_matcher.initialize(&m_g, &m_av_manager, t_n_threads);
            // sc.initialize(&g);

            // partitioning
            m_kaffpa_partitioner.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, m_hierarchy, m_distance, t_imbalance, t_n_threads);

            // refinement
            // m_i_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, m_hierarchy, m_distance, m_lmax);
            m_lp_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, &m_qgraph, m_hierarchy, m_distance, m_lmax, t_n_threads);
            // qgr.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o);

            auto ep_io = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            internal_solve();

#if STATISTICCOLLECTOR
            weight_t qap = get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle);
            std::vector<weight_t> pweights;
            for (partition_t id = 0; id < m_k; ++id) { pweights.push_back(m_p_manager.get_bweight(id)); }
            m_stat_collect.set_final(qap, pweights, m_lmax);
#endif
            m_stat_collect.finalize();

            std::cout << m_stat_collect.to_JSON() << std::endl;

            std::vector<partition_t> p(m_g.get_n());
            for (vertex_t u = 0; u < m_g.get_n(); ++u) { p[u] = m_p_manager[u]; }

            return p;
        }

    private:
        void internal_solve() {
            s32 level = 0;

            while (!(m_av_manager.get_n_active() <= m_threshold || m_av_manager.get_n_active() <= m_k * 64)) {
                matching(level);
                coarsening(level);
                level += 1;
            }

            partition();

            while (level > 0) {
                level -= 1;
                uncoarsening(level);
                refinement(level);
            }
        }

        void partition() {
            auto sp_partition = std::chrono::high_resolution_clock::now();

            m_kaffpa_partitioner.partition();

            auto ep_partition = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));

            // HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_partition_stats(get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle), m_p_manager.get_bweights(), m_lmax);
#endif
        }

        void matching(s32 level) {
            auto sp_match = std::chrono::high_resolution_clock::now();

            m_matches.emplace_back();
            m_matches.back().reserve(m_av_manager.get_n_active() / 2);

            m_he_matcher.match(m_matches.back());

            auto ep_match = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);

#if STATISTICCOLLECTOR
            m_stat_collect.set_matching_stats(level, m_matches.back().size());
#endif
        }

        void coarsening(s32 level) {
            auto sp_coarse = std::chrono::high_resolution_clock::now();

// #pragma omp parallel for default(none) num_threads(m_n_threads)
            for (auto &e: m_matches.back()) {
                m_g.contract(e.u, e.v);
                m_av_manager.contract(e.u, e.v);
            }

            auto ep_coarse = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);

            // HEAVYASSERT(assert_state_pre_partitioning(m_g, m_av_manager));
#if STATISTICCOLLECTOR
            m_stat_collect.set_coarsening_stats(m_av_manager.get_n_active(), level);
#endif
        }

        void uncoarsening(s32 level) {
            auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            for (size_t i = 0; i < m_matches.back().size(); ++i) {
                u64 idx = m_matches.back().size() - 1 - i;
                vertex_t u = m_matches.back()[idx].u;
                vertex_t v = m_matches.back()[idx].v;

                m_g.uncontract(u, v);
                m_av_manager.uncontract(u, v);
                m_p_manager.uncontract(u, v);
                m_bv_manager.uncontract(u, v);
            }
            m_matches.pop_back();

            auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_uncoarsening_time(get_seconds(sp_uncoarse, ep_uncoarse), level);

            // HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_uncoarsening_stats(level, m_av_manager.get_n_active());
#endif
        }

        void refinement(s32 level) {
            auto sp_refinement = std::chrono::high_resolution_clock::now();

            // m_i_refine.refine();
            m_lp_refine.refine();
            // qgr.refine(pm);

            auto ep_refinement = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), level);

            HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_refinement_stats(level, get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle));
#endif
        }
    };
}
#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_SOLVER_H
