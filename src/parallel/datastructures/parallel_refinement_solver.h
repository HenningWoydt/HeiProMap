#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H

#include "../utility/parallel_utils.h"

#include "parallel_graph.h"
#include "../../serial/datastructures/active_vertex_manager.h"
#include "../../serial/datastructures/partition_manager.h"
#include "../../serial/datastructures/boundary_vertex_manager.h"
#include "../../serial/refinement/label_propagation_refinement.h"
#include "../../serial/refinement/identity_refinement.h"
#include "../../serial/datastructures/statistic_collector.h"
#include "parallel_distance_oracle.h"
#include "parallel_partition_manager.h"
#include "parallel_active_vertex_manager.h"
#include "parallel_boundary_vertex_manager.h"
#include "parallel_static_csr_graph.h"
#include "../refinement/parallel_label_propagation_refinement.h"
#include "parallel_quotient_graph.h"
#include "../refinement/parallel_quotient_graph_refinement.h"


namespace HeiProMap {

    class ParallelRefinementSolver {
    private:
        // main structures
        std::string m_mapping_in;
        ParallelStaticCSRGraph m_g;
        ParallelActiveVertexManager< typeof(m_g) > m_av_manager;
        ParallelPartitionManager<typeof(m_g), typeof(m_av_manager)> m_p_manager;
        ParallelBoundaryVertexManager<typeof(m_g), typeof(m_av_manager), typeof(m_p_manager)> m_bv_manager;
        ParallelDistanceOracle m_d_oracle;
        ParallelQuotientGraph<typeof(m_g), typeof(m_av_manager), typeof(m_p_manager), typeof(m_d_oracle)> m_qgraph;

        // distance
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k;

        // mult threading
        u64 m_n_threads;

        // balance
        weight_t m_lmax = 0;

        // refinement
        ParallelLabelPropagationRefinement<typeof(m_g), typeof(m_av_manager), typeof(m_bv_manager), typeof(m_p_manager), typeof(m_d_oracle), typeof(m_qgraph)> m_lp_refine;
        // IdentityRefinement m_i_refine;
        ParallelQuotientGraphRefinement<typeof(m_g), typeof(m_av_manager), typeof(m_bv_manager), typeof(m_p_manager), typeof(m_d_oracle), typeof(m_qgraph)> m_qg_refine;

        // statistics
        StatisticCollector m_stat_collect;

    public:
        ParallelRefinementSolver(const std::string &t_graph_in,
                                 const std::string &t_mapping_in,
                                 std::vector<partition_t> &t_hierarchy,
                                 std::vector<weight_t> &t_distance,
                                 f64 t_imbalance,
                                 u64 n_threads) {
            auto sp_graph_io = std::chrono::high_resolution_clock::now();
            m_g.initialize(t_graph_in, n_threads);
            auto ep_graph_io = std::chrono::high_resolution_clock::now();

            auto sp_io = std::chrono::high_resolution_clock::now();

            // distance
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_k = prod<partition_t>(m_hierarchy);

            // manager
            m_mapping_in = t_mapping_in;
            m_av_manager.initialize(&m_g, n_threads);
            m_p_manager.initialize(&m_g, &m_av_manager, m_k, n_threads);
            m_bv_manager.initialize(&m_g, &m_av_manager, &m_p_manager, m_k, n_threads);
            m_d_oracle.initialize(m_hierarchy, m_distance, n_threads);
            m_qgraph.initialize(&m_g, &m_p_manager, &m_d_oracle, m_k, n_threads);
            HEAVYASSERT(assert_state_pre_partitioning(m_g, m_av_manager));

            // mult threading
            m_n_threads = n_threads;

            // balance
            m_lmax = ceil((1.0 + t_imbalance) * ((f64) m_g.get_weight() / (f64) m_k));

            // refinement
            // m_i_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, m_hierarchy, m_distance, m_lmax);
            m_lp_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, &m_qgraph, m_hierarchy, m_distance, m_lmax, n_threads);
            m_qg_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, &m_qgraph, m_hierarchy, m_distance, m_lmax, n_threads);

            auto ep_io = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            partition();
            refine();

#if STATISTICCOLLECTOR
            weight_t qap = get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle);
            std::vector<weight_t> pweights;
            for (partition_t id = 0; id < m_k; ++id) { pweights.push_back(m_p_manager.get_bweight(id)); }
            m_stat_collect.set_final(qap, pweights, m_lmax);
#endif
            m_stat_collect.finalize();

            std::cout << m_stat_collect.to_JSON() << std::endl;

            /*
            std::vector<partition_t> p(m_g.get_n());
#pragma omp parallel for default(none) shared(p) num_threads(m_n_threads)
            for (vertex_t u = 0; u < m_g.get_n(); ++u) { p[u] = m_p_manager[u]; }
             */

            return {}; //p;
        }

    private:
        void partition() {
            auto sp_partition = std::chrono::high_resolution_clock::now();

            // read in the vector
            std::vector<partition_t> p(m_g.get_n());
            if (m_n_threads == 1) {
                read_partition(m_mapping_in, p);
            } else {
                parallel_read_partition(m_mapping_in, p, m_n_threads);
            }

#pragma omp parallel for default(none) shared(p) num_threads(m_n_threads)
            for(vertex_t u = 0; u < m_g.get_n(); ++u){
                m_p_manager.set(u, p[u]);
            }

            // then initialize boundary vertices
#pragma omp parallel default(none) num_threads(m_n_threads)
            {
                size_t thread_id = omp_get_thread_num();
                partition_t base_range = floor((f64) m_k / (f64) m_n_threads);
                partition_t rem = m_k % m_n_threads;

                partition_t start_b;
                partition_t end_b;
                if (thread_id < rem) {
                    start_b = thread_id * (base_range + 1);
                    end_b = start_b + base_range + 1;
                } else {
                    start_b = rem * (base_range + 1) + (thread_id - rem) * base_range;
                    end_b = start_b + base_range;
                }

                for (vertex_t u = 0; u < m_g.get_n(); ++u) {
                    if (start_b <= m_p_manager[u] && m_p_manager[u] < end_b) {
                        m_bv_manager.insert(u, m_p_manager[u]);
                    }
                }
            }
            
            // then initialize quotient graph
            for(vertex_t u = 0; u < m_g.get_n(); ++u){
                for(size_t i = 0; i < m_g.size(u); ++i){
                    vertex_t v = m_g.neighbor(u, i);
                    weight_t w = m_g.get_weight(u, i);

                    partition_t u_id = m_p_manager[u];
                    partition_t v_id = m_p_manager[v];

                    m_qgraph.add_edge(u_id, v_id, w);
                }
            }

            auto ep_partition = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_partition_stats(get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle), m_p_manager.get_bweights(), m_lmax);
#endif
        }

        void refine() {
            auto sp_refinement = std::chrono::high_resolution_clock::now();

            for(size_t i = 0; i < 1; ++i) {
                // m_i_refine.refine();
                m_lp_refine.refine();
                HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
                m_qg_refine.refine();
                HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
            }

            auto ep_refinement = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), 0);

            HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_refinement_stats(0, get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle));
#endif
        }
    };

}

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H
