#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H

namespace HeiProMap {

    class ParallelRefinementSolver {
    private:
        // main structures
        Graph m_g;
        ActiveVertexManager m_av_manager;
        PartitionManager m_p_manager;
        BoundaryVertexManager m_bv_manager;

        // distance
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k;
        DistanceOracle m_d_oracle;

        // balance
        weight_t m_lmax = 0;

        // multilevel
        vertex_t m_threshold;

        // matching
        // GreedyEdgeMatcher gem;
        HeavyEdgeMatcher m_he_matcher;
        // SimpleClustering sc;

        std::vector<std::vector<Edge>> m_matches;

        // partitioning
        KaffpaPartitioner m_kaffpa_partitioner;

        // refinement
        LabelPropagationRefinement m_lp_refine;
        IdentityRefinement m_i_refine;
        // QuotientGraphRefinement qgr;

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
            m_g.initialize(t_graph_in, 1);
            auto ep_graph_io = std::chrono::high_resolution_clock::now();

            auto sp_io = std::chrono::high_resolution_clock::now();

            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_k = prod<partition_t>(m_hierarchy);

            // manager
            m_av_manager.initialize(&m_g);
            m_p_manager.initialize(&m_g, &m_av_manager, m_k);
            m_bv_manager.initialize(&m_g, &m_av_manager, &m_p_manager, m_k);
            HEAVYASSERT(assert_state_pre_partitioning(m_g, m_av_manager));

            // distance
            m_d_oracle.initialize(m_hierarchy, m_distance);

            // balance
            m_lmax = ceil((1.0 + t_imbalance) * ((f64) m_g.get_weight() / (f64) m_k));

            // multilevel
            m_threshold = m_g.get_n() / 500;

            // matching
            // gem.initialize(&g);
            m_he_matcher.initialize(&m_g, &m_av_manager);
            // sc.initialize(&g);

            // partitioning
            m_kaffpa_partitioner.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, m_hierarchy, m_distance, t_imbalance);

            // refinement
            m_i_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, m_hierarchy, m_distance, m_lmax);
            m_lp_refine.initialize(&m_g, &m_av_manager, &m_bv_manager, &m_p_manager, &m_d_oracle, m_hierarchy, m_distance, m_lmax);
            // qgr.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o);

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

            std::vector<partition_t> p(m_g.get_n());
            for (vertex_t u = 0; u < m_g.get_n(); ++u) { p[u] = m_p_manager[u]; }

            return p;
        }

    private:
        void partition() {
            auto sp_partition = std::chrono::high_resolution_clock::now();

            m_kaffpa_partitioner.partition();

            auto ep_partition = std::chrono::high_resolution_clock::now();
            m_stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(m_g, m_av_manager, m_p_manager, m_bv_manager, m_k));
#if STATISTICCOLLECTOR
            m_stat_collect.set_partition_stats(get_qap(m_g, m_av_manager, m_p_manager, m_d_oracle), m_p_manager.get_bweights(), m_lmax);
#endif
        }

        void refine() {
            auto sp_refinement = std::chrono::high_resolution_clock::now();

            // m_i_refine.refine();
            m_lp_refine.refine();
            // qgr.refine(pm);

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
