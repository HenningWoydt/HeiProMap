#ifndef HEIDELBERGPROCESSMAPPING_SOLVER_H
#define HEIDELBERGPROCESSMAPPING_SOLVER_H

#include "active_vertex_manager.h"
#include "boundary_vertex_manager.h"
#include "csr_graph.h"
#include "graph.h"
#include "partition_manager.h"
#include "statistic_collector.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/simple_clustering.h"
#include "../coarsening/simple_edge_matcher.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/identity_refinement.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class Solver {
        CSRGraph g;
        ActiveVertexManager<typeof(g)> av_manager;
        PartitionManager<typeof(g), typeof(av_manager)> p_manager;
        BoundaryVertexManager<typeof(g), typeof(av_manager), typeof(p_manager)> bv_manager;

        // distance
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k;
        DistanceOracle d_oracle;

        // balance
        weight_t lmax = 0;

        // multilevel
        vertex_t threshold;

        // matching
        GreedyEdgeMatcher<typeof(g), typeof(av_manager)> ge_matcher;
        HeavyEdgeMatcher<typeof(g), typeof(av_manager)> he_matcher;
        // SimpleClustering sc;

        std::vector<std::vector<EdgeUV>> matches;

        // partitioning
        KaffpaPartitioner<typeof(g), typeof(av_manager), typeof(bv_manager), typeof(p_manager)> kaffpa_partitioner;

        // refinement
        LabelPropagationRefinement<typeof(g), typeof(av_manager), typeof(bv_manager), typeof(p_manager), typeof(d_oracle)> lp_refine;
        IdentityRefinement<typeof(g), typeof(av_manager), typeof(bv_manager), typeof(p_manager), typeof(d_oracle)> i_refine;
        // QuotientGraphRefinement qgr;

        // statistics
        StatisticCollector stat_collect;

    public:
        Solver(const std::string& t_graph_in,
               const std::vector<partition_t>& t_hierarchy,
               const std::vector<weight_t>& t_distance,
               const f64 t_imbalance) {
            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            g.initialize(t_graph_in);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();

            const auto sp_io = std::chrono::high_resolution_clock::now();

            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);

            // manager
            av_manager.initialize(&g);
            p_manager.initialize(&g, &av_manager, k);
            bv_manager.initialize(&g, &av_manager, &p_manager, k);
            HEAVYASSERT(assert_state_pre_partitioning(g, av_manager));

            // distance
            d_oracle.initialize(hierarchy, distance);

            // balance
            lmax = ceil((1.0 + t_imbalance) * ((f64)g.get_weight() / (f64)k));

            // multilevel
            threshold = g.get_n() / 500;

            // matching
            ge_matcher.initialize(&g, &av_manager);
            he_matcher.initialize(&g, &av_manager);
            // sc.initialize(&g);

            // partitioning
            kaffpa_partitioner.initialize(&g, &av_manager, &bv_manager, &p_manager, hierarchy, distance, t_imbalance);

            // refinement
            i_refine.initialize(&g, &av_manager, &bv_manager, &p_manager, &d_oracle, hierarchy, distance, lmax);
            lp_refine.initialize(&g, &av_manager, &bv_manager, &p_manager, &d_oracle, hierarchy, distance, lmax);
            // qgr.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o);

            const auto ep_io = std::chrono::high_resolution_clock::now();
            stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            internal_solve();

#if STATISTICCOLLECTOR
            weight_t qap = get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle);
            std::vector<weight_t> pweights;
            for (partition_t id = 0; id < m_k; ++id) { pweights.push_back(m_p_manager.get_bweight(id)); }
            m_stat_collect.set_final(qap, pweights, m_lmax);
#endif
            stat_collect.finalize();

            std::cout << stat_collect.to_JSON() << std::endl;

            std::vector<partition_t> p(g.get_n());
            for (vertex_t u = 0; u < g.get_n(); ++u) { p[u] = p_manager[u]; }

            return p;
        }

    private:
        void internal_solve() {
            s32 level = 0;

            while (!(av_manager.get_n_active() <= threshold || av_manager.get_n_active() <= k * 64)) {
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
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            kaffpa_partitioner.partition();

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(g, av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_partition_stats(get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle), m_p_manager.get_bweights(), m_lmax);
#endif
        }

        void matching(s32 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().reserve(av_manager.get_n_active() / 2);

            // ge_matcher.match(matches.back());
            he_matcher.match(matches.back());

            const auto ep_match = std::chrono::high_resolution_clock::now();
            stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);

#if STATISTICCOLLECTOR
            m_stat_collect.set_matching_stats(level, m_matches.back().size());
#endif
        }

        void coarsening(s32 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            for (const auto& e : matches.back()) {
                g.contract(e.u, e.v);
                av_manager.contract(e.u, e.v);
            }

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);

            HEAVYASSERT(assert_state_pre_partitioning(g, av_manager));
#if STATISTICCOLLECTOR
            m_stat_collect.set_coarsening_stats(m_av_manager.get_n_active(), level);
#endif
        }

        void uncoarsening(s32 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            for (size_t i = 0; i < matches.back().size(); ++i) {
                const u64 idx    = matches.back().size() - 1 - i;
                const vertex_t u = matches.back()[idx].u;
                const vertex_t v = matches.back()[idx].v;

                g.uncontract(u, v);
                av_manager.uncontract(u, v);
                p_manager.uncontract(u, v);
                bv_manager.uncontract(u, v);
            }
            matches.pop_back();

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            stat_collect.set_uncoarsening_time(get_seconds(sp_uncoarse, ep_uncoarse), level);

            HEAVYASSERT(assert_state_after_partitioning(g, av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_uncoarsening_stats(level, m_av_manager.get_n_active());
#endif
        }

        void refinement(s32 level) {
            const auto sp_refinement = std::chrono::high_resolution_clock::now();

            i_refine.refine();
            lp_refine.refine();
            // qgr.refine(pm);

            const auto ep_refinement = std::chrono::high_resolution_clock::now();
            stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), level);

            HEAVYASSERT(assert_state_after_partitioning(g, av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_refinement_stats(level, get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle));
#endif
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_SOLVER_H
