#ifndef SERIALPROCESSMAPPING_STATISTIC_COLLECTOR_H
#define SERIALPROCESSMAPPING_STATISTIC_COLLECTOR_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/JSON_utils.h"

#ifndef STATISTICCOLLECTOR
#define STATISTICCOLLECTOR true
#endif

namespace HeiProMap {

    class StatisticCollector {

    private:
        u64 max_level = 0;
        f64 total_time = 0.0;

        // io
        f64 total_graph_io_time = 0.0;
        f64 total_io_time = 0.0;

        // matching
        f64 total_matching_time = 0.0;
#if STATISTICCOLLECTOR
        std::vector<f64> matching_time;
        std::vector<vertex_t> matching_size;
        vertex_t total_matching_size = 0;
#endif

        // coarsening
#if STATISTICCOLLECTOR
        std::vector<f64> coarsening_time;
        std::vector<vertex_t> coarsening_end_size;
#endif
        f64 total_coarsening_time = 0.0;

        // uncoarsening
#if STATISTICCOLLECTOR
        std::vector<f64> uncoarsening_time;
        std::vector<vertex_t> uncoarsening_end_size;
#endif
        f64 total_uncoarsening_time = 0.0;

        // refinement
        f64 total_refinement_time = 0.0;
#if STATISTICCOLLECTOR
        std::vector<f64> refinement_time;
        std::vector<u64> refinement_end_qap;
#endif

        // partition
        f64 totaL_partition_time = 0.0;
#if STATISTICCOLLECTOR
        u64 partition_end_qap = 0;
        weight_t partition_end_max_pweight = 0.0;
        f64 partition_end_avg_pweight = 0.0;
        weight_t partition_end_min_pweight = 0.0;
        std::vector<weight_t> partition_end_pweights;
        bool partition_end_balanced = false;
#endif

        // final information
#if STATISTICCOLLECTOR
        u64 final_qap = 0;
        weight_t final_max_pweight = 0.0;
        f64 final_avg_pweight = 0.0;
        weight_t final_min_pweight = 0.0;
        std::vector<weight_t> final_pweights;
        u64 final_lmax = 0;
        bool final_balanced = false;
#endif


    public:
        StatisticCollector() {
            resize(1);
            max_level = 0;
        }


        inline void finalize() {
            total_time = total_graph_io_time + total_io_time + total_matching_time + total_coarsening_time + totaL_partition_time + total_uncoarsening_time + total_refinement_time;
        }

        inline void set_io(f64 graph_time, f64 time) {
            total_graph_io_time = graph_time;
            total_io_time = time;
        }

        inline void set_matching_time(f64 time, u64 level) {
            total_matching_time += time;
#if STATISTICCOLLECTOR
            resize(level);
            matching_time[level] = time;
#endif
        }

        inline void set_matching_stats(u64 level, vertex_t size) {
#if STATISTICCOLLECTOR
            matching_size[level] = size;
            total_matching_size += size;
#endif
        }


        inline void set_coarsening_time(f64 time, u64 level) {
            total_coarsening_time += time;
#if STATISTICCOLLECTOR
            resize(level);
            coarsening_time[level] = time;
#endif
        }

        inline void set_coarsening_stats(vertex_t size, u64 level) {
#if STATISTICCOLLECTOR
            coarsening_end_size[level] = size;
#endif
        }

        inline void set_uncoarsening_time(f64 time, u64 level) {
            total_uncoarsening_time += time;
#if STATISTICCOLLECTOR
            resize(level);
            uncoarsening_time[level] = time;
#endif
        }

        inline void set_uncoarsening_stats(u64 level, vertex_t size) {
#if STATISTICCOLLECTOR
            uncoarsening_end_size[level] = size;
#endif
        }

        inline void set_refinement_time(f64 time, u64 level) {
            total_refinement_time += time;
#if STATISTICCOLLECTOR
            resize(level);
            refinement_time[level] = time;
#endif
        }

        inline void set_refinement_stats(u64 level, u64 objective) {
#if STATISTICCOLLECTOR
            refinement_end_qap[level] = objective;
#endif
        }

        inline void set_partition_time(f64 time) { totaL_partition_time = time; }

        inline void set_partition_stats(u64 end_objective, const std::vector<weight_t> &pweights, weight_t lmax) {
#if STATISTICCOLLECTOR
            partition_end_qap = end_objective;
            partition_end_pweights = pweights;
            partition_end_max_pweight = max(pweights);
            partition_end_avg_pweight = (f64) avg(pweights);
            partition_end_min_pweight = min(pweights);
            partition_end_balanced = partition_end_max_pweight <= lmax;
#endif
        }

        inline void set_final(u64 end_objective, const std::vector<weight_t> &pweights, weight_t lmax) {
#if STATISTICCOLLECTOR
            final_qap = end_objective;
            final_pweights = pweights;
            final_max_pweight = max(pweights);
            final_avg_pweight = (f64) avg(pweights);
            final_min_pweight = min(pweights);
            final_lmax = lmax;
            final_balanced = final_max_pweight <= lmax;
#endif
        }

        std::string to_JSON() {
            std::string s = "{\n";

            s += to_JSON_MACRO(total_time);
            s += to_JSON_MACRO(total_graph_io_time);
            s += to_JSON_MACRO(total_io_time);
            s += to_JSON_MACRO(total_matching_time);
            s += to_JSON_MACRO(total_coarsening_time);
            s += to_JSON_MACRO(totaL_partition_time);
            s += to_JSON_MACRO(total_uncoarsening_time);
            s += to_JSON_MACRO(total_refinement_time);
#if STATISTICCOLLECTOR
            s += to_JSON_MACRO(max_level);
            s += to_JSON_MACRO(matching_time);
            s += to_JSON_MACRO(matching_size);
            s += to_JSON_MACRO(coarsening_time);
            s += to_JSON_MACRO(coarsening_end_size);
            s += to_JSON_MACRO(uncoarsening_time);
            s += to_JSON_MACRO(uncoarsening_end_size);
            s += to_JSON_MACRO(refinement_time);
            s += to_JSON_MACRO(refinement_end_qap);
            s += to_JSON_MACRO(partition_end_qap);
            s += to_JSON_MACRO(partition_end_max_pweight);
            s += to_JSON_MACRO(partition_end_avg_pweight);
            s += to_JSON_MACRO(partition_end_min_pweight);
            s += to_JSON_MACRO(partition_end_pweights);
            s += to_JSON_MACRO(partition_end_balanced);
            s += to_JSON_MACRO(final_qap);
            s += to_JSON_MACRO(final_max_pweight);
            s += to_JSON_MACRO(final_avg_pweight);
            s += to_JSON_MACRO(final_min_pweight);
            s += to_JSON_MACRO(final_pweights);
            s += to_JSON_MACRO(final_lmax);
            s += to_JSON_MACRO(final_balanced);
#endif

            s.pop_back();
            s.pop_back();
            s += "\n}";
            return s;
        }

    private:
        void resize(u64 new_max_level) {
#if STATISTICCOLLECTOR
            if (new_max_level > max_level) {
                // matching
                matching_time.resize(new_max_level + 1, 0.0);
                matching_size.resize(new_max_level + 1, 0);

                // coarsening
                coarsening_time.resize(new_max_level + 1, 0.0);
                coarsening_end_size.resize(new_max_level + 1, 0);

                // coarsening
                uncoarsening_time.resize(new_max_level + 1, 0.0);
                uncoarsening_end_size.resize(new_max_level + 1, 0);

                // refinement
                refinement_time.resize(new_max_level + 1, 0.0);
                refinement_end_qap.resize(new_max_level + 1, 0);

                max_level = new_max_level;
            }
#endif
        }

    };

}

#endif //SERIALPROCESSMAPPING_STATISTIC_COLLECTOR_H
