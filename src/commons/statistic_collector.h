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

#ifndef HEIPROMAP_STATISTIC_COLLECTOR_H
#define HEIPROMAP_STATISTIC_COLLECTOR_H

#ifndef COLLECT_METRICS
#define COLLECT_METRICS true
#endif

#if COLLECT_METRICS
#define METRICS(x) x
#else
#define METRICS(x)
#endif

#if COLLECT_METRICS
#define METRICS_TIME(x) const auto x = std::chrono::high_resolution_clock::now();
#else
#define METRICS_TIME(x)
#endif

#include "JSON_utils.h"

namespace HeiProMap {
    class StatisticCollector {
#if COLLECT_METRICS
        u64 max_level  = 0;
        f64 total_time = 0.0;

        // io
        f64 total_graph_io_time = 0.0;
        f64 total_io_time       = 0.0;

        // matching
        f64 total_matching_time = 0.0;
        std::vector<f64> matching_time      = {0};
        std::vector<vertex_t> matching_size = {0};
        vertex_t total_matching_size        = 0;
        JSONString matching_method_stats;

        // coarsening
        f64 total_coarsening_time = 0.0;
        std::vector<f64> coarsening_time          = {0.0};
        std::vector<vertex_t> coarsening_end_size = {0};

        // uncoarsening
        f64 total_uncoarsening_time = 0.0;
        std::vector<f64> uncoarsening_time          = {0.0};
        std::vector<vertex_t> uncoarsening_end_size = {0};

        // refinement
        f64 total_refinement_time = 0.0;
        std::map<std::string, JSONString> refinement_method_stats;
        std::vector<f64> refinement_time    = {0.0};
        std::vector<u64> refinement_end_qap = {0};

        // partition
        f64 total_partition_time = 0.0;
        u64 partition_end_qap              = 0;
        weight_t partition_end_max_pweight = 0.0;
        f64 partition_end_avg_pweight      = 0.0;
        weight_t partition_end_min_pweight = 0.0;
        std::vector<weight_t> partition_end_pweights;
        bool partition_end_balanced = false;

        // final information
        u64 final_qap              = 0;
        weight_t final_max_pweight = 0.0;
        f64 final_avg_pweight      = 0.0;
        weight_t final_min_pweight = 0.0;
        std::vector<weight_t> final_pweights;
        u64 final_lmax      = 0;
        bool final_balanced = false;

    public:
        StatisticCollector() {
            max_level = 0;
        }


        void finalize() {
            total_time = total_graph_io_time + total_io_time + total_matching_time + total_coarsening_time + total_partition_time + total_uncoarsening_time + total_refinement_time;
        }

        void set_io(const f64 graph_time, const f64 time) {
            total_graph_io_time = graph_time;
            total_io_time       = time;
        }

        void set_matching_time(const f64 time, const u64 level) {
            total_matching_time += time;
            resize(level);
            matching_time[level] = time;
        }

        void set_matching_stats(const u64 level, const vertex_t size) {
            matching_size[level] = size;
            total_matching_size += size;
        }


        void set_coarsening_time(const f64 time, const u64 level) {
            total_coarsening_time += time;
            resize(level);
            coarsening_time[level] = time;
        }

        void set_coarsening_stats(const vertex_t size, const u64 level) {
            coarsening_end_size[level] = size;
        }

        void set_uncoarsening_time(const f64 time, const u64 level) {
            total_uncoarsening_time += time;
            resize(level);
            uncoarsening_time[level] = time;
        }

        void set_uncoarsening_stats(const u64 level, const vertex_t size) {
            uncoarsening_end_size[level] = size;
        }

        void set_refinement_time(const f64 time, const u64 level) {
            total_refinement_time += time;
            resize(level);
            refinement_time[level] = time;
        }

        void set_refinement_stats(const u64 level, const u64 objective) {
            refinement_end_qap[level] = objective;
        }

        void add_refinement_method_stats(const std::string& method, const JSONString& stats) {
            refinement_method_stats.emplace(method, stats);
        }

        void add_matching_method_stats(const JSONString& stats) {
            matching_method_stats = stats;
        }

        void set_partition_time(const f64 time) { total_partition_time = time; }

        void set_partition_stats(const u64 end_objective, const std::vector<weight_t>& pweights, const weight_t lmax) {
            partition_end_qap         = end_objective;
            partition_end_pweights    = pweights;
            partition_end_max_pweight = max(pweights);
            partition_end_avg_pweight = (f64)avg(pweights);
            partition_end_min_pweight = min(pweights);
            partition_end_balanced    = partition_end_max_pweight <= lmax;
        }

        void set_final(const u64 end_objective, const std::vector<weight_t>& pweights, const weight_t lmax) {
            final_qap         = end_objective;
            final_pweights    = pweights;
            final_max_pweight = max(pweights);
            final_avg_pweight = (f64)avg(pweights);
            final_min_pweight = min(pweights);
            final_lmax        = lmax;
            final_balanced    = final_max_pweight <= lmax;
        }

        std::string to_JSON() const {
            std::string s = "{ \n";

            s += to_JSON_MACRO(total_time);
            s += to_JSON_MACRO(total_graph_io_time);
            s += to_JSON_MACRO(total_io_time);
            s += to_JSON_MACRO(total_matching_time);
            s += to_JSON_MACRO(total_coarsening_time);
            s += to_JSON_MACRO(total_partition_time);
            s += to_JSON_MACRO(total_uncoarsening_time);
            s += to_JSON_MACRO(total_refinement_time);
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

            s += to_JSON_MACRO(matching_method_stats);
            s += to_JSON_MACRO(refinement_method_stats);

            s.pop_back();
            s.pop_back();
            s += "\n}";
            return s;
        }

    private:
        void resize(const u64 new_max_level) {
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
        }
#endif
    };
}

#endif //HEIPROMAP_STATISTIC_COLLECTOR_H
