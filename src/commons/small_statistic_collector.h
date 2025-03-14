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

#ifndef HEIPROMAP_SMALL_STATISTIC_COLLECTOR_H
#define HEIPROMAP_SMALL_STATISTIC_COLLECTOR_H

#ifndef COLLECT_SMALL_METRICS
#define COLLECT_SMALL_METRICS true
#endif

#if COLLECT_SMALL_METRICS
#define SMALL_METRICS(x) x
#else
#define SMALL_METRICS(x)
#endif

#include <string>
#include <map>
#include <iostream>
#include <iomanip>  // for std::setw, std::setprecision
#include <map>
#include <string>
#include <utility> // for std::pair

#include "../definitions.h"

namespace HeiProMap {

    class SmallStatisticCollector {
        std::map<std::string, f64>                 method_time;
        std::map<std::string, std::pair<f64, s64>> refinement_stats;

    public:
        void add(std::string method, f64 time) {
            method_time[method] += time;
        }

        void add_refinement(std::string method, f64 time, s64 qap_delta) {
            auto &ref_stats = refinement_stats[method];
            ref_stats.first += time;
            ref_stats.second += qap_delta;
        }

        void print() {
            print_method_times();
            print_refinement_stats();
        }

    private:
        // Helper for printing the "Method Times" table
        void print_method_times() {
            // 1) Compute total
            double          totalTime = 0.0;
            for (const auto &kv: method_time) {
                totalTime += kv.second;
            }

            // 2) Prepare rows (header + data rows + total row)
            // Each row is a vector of strings: {MethodName, Time(s), Time(%)}
            std::vector<std::vector<std::string>> rows;

            // Header row
            rows.push_back({"Method", "Time(s)", "Time(%)"});

            // Data rows
            for (const auto &kv: method_time) {
                const std::string  &method = kv.first;
                double             time    = kv.second;
                // Build strings for time and percentage
                std::ostringstream time_str;
                time_str << std::fixed << std::setprecision(4) << time;

                double             pct = (totalTime > 0.0) ? (time / totalTime * 100.0) : 0.0;
                std::ostringstream pct_str;
                pct_str << std::fixed << std::setprecision(2) << pct << "%";

                rows.push_back({method, time_str.str(), pct_str.str()});
            }

            // Total row
            {
                std::ostringstream total_time_str;
                total_time_str << std::fixed << std::setprecision(4) << totalTime;
                rows.push_back({"Total", total_time_str.str(), "100.00%"});
            }

            // 3) Compute column widths
            // We have 3 columns: method, time(s), time(%)
            const int           numCols = 3;
            std::vector<size_t> colWidths(numCols, 0);
            for (const auto     &row: rows) {
                for (int c = 0; c < numCols; ++c) {
                    // Update max width by comparing current stored width to this entry's length
                    if (row[c].size() > colWidths[c]) {
                        colWidths[c] = row[c].size();
                    }
                }
            }

            // 4) Print the table
            std::cout << "\n=== Method Times ===\n";

            // Print header row
            printRow(rows[0], colWidths, /*isHeader=*/true);

            // Print a separator line
            printSeparator(colWidths);

            // Print data rows
            for (size_t i = 1; i < rows.size(); i++) {
                printRow(rows[i], colWidths);
            }
        }

        // Helper for printing the "Refinement Stats" table
        void print_refinement_stats() {
            // 1) Calculate totals
            double          totalRefTime  = 0.0;
            s64             totalQAPDelta = 0;
            for (const auto &kv: refinement_stats) {
                totalRefTime += kv.second.first;
                totalQAPDelta += kv.second.second;
            }

            // 2) Prepare rows
            // Now we have: { "Refinement", "Time(s)", "Time(%)", "QAP", "QAP(%)", "QAP/Time" }
            std::vector<std::vector<std::string>> rows;
            rows.push_back({"Refinement", "Time(s)", "Time(%)", "QAP", "QAP(%)", "QAP/Time"}); // header

            // Data rows
            for (const auto &kv: refinement_stats) {
                const auto &method  = kv.first;
                double     time     = kv.second.first;
                s64        qapDelta = kv.second.second;

                // Time(s)
                std::ostringstream time_str;
                time_str << std::fixed << std::setprecision(4) << time;

                // Time(%)
                double             pct = (totalRefTime > 0.0) ? (time / totalRefTime * 100.0) : 0.0;
                std::ostringstream pct_str;
                pct_str << std::fixed << std::setprecision(2) << pct << "%";

                // QAP
                std::ostringstream qap_str;
                qap_str << qapDelta;

                // QAP(%)
                double             qapPct = (totalQAPDelta != 0)
                                            ? (static_cast<double>(qapDelta) / totalQAPDelta * 100.0)
                                            : 0.0;
                std::ostringstream qap_pct_str;
                qap_pct_str << std::fixed << std::setprecision(2) << qapPct << "%";

                // QAP / Time
                double eff = 0.0;
                if (time > 0.0) {
                    eff = qapDelta / time;
                }
                std::ostringstream eff_str;
                // You can choose the precision to display here
                eff_str << std::fixed << std::setprecision(2) << eff;

                rows.push_back({
                                       method,
                                       time_str.str(),
                                       pct_str.str(),
                                       qap_str.str(),
                                       qap_pct_str.str(),
                                       eff_str.str()
                               });
            }

            // Total row
            {
                std::ostringstream time_str;
                time_str << std::fixed << std::setprecision(4) << totalRefTime;

                std::ostringstream qap_str;
                qap_str << totalQAPDelta;

                // For the total row, QAP/Time might be somewhat ambiguous,
                // but you can define it if you want:
                double total_eff = 0.0;
                if (totalRefTime > 0.0) {
                    total_eff = static_cast<double>(totalQAPDelta) / totalRefTime;
                }
                std::ostringstream total_eff_str;
                total_eff_str << std::fixed << std::setprecision(2) << total_eff;

                rows.push_back({
                                       "Total",
                                       time_str.str(),
                                       "100.00%",
                                       qap_str.str(),
                                       "100.00%",
                                       total_eff_str.str()
                               });
            }

            // 3) Determine column widths
            auto colWidths = computeColumnWidths(rows);

            // 4) Print the table
            std::cout << "\n=== Refinement Stats ===\n";
            printRow(rows[0], colWidths, true);    // header
            printSeparator(colWidths);
            for (size_t i = 1; i < rows.size(); i++) {
                printRow(rows[i], colWidths, false);
            }

            std::cout << std::endl;
        }

        // Computes the maximum width needed for each column across all rows
        std::vector<size_t> computeColumnWidths(const std::vector<std::vector<std::string>> &rows) {
            if (rows.empty()) return {};

            size_t              numCols = rows[0].size();
            std::vector<size_t> colWidths(numCols, 0);

            for (auto &row: rows) {
                for (size_t c = 0; c < numCols; ++c) {
                    if (row[c].size() > colWidths[c]) {
                        colWidths[c] = row[c].size();
                    }
                }
            }
            return colWidths;
        }

        // Utility to print a single row of data, given the precomputed column widths
        void printRow(const std::vector<std::string> &row, const std::vector<size_t> &colWidths, bool isHeader = false) {
            // We'll left-align the first column (name) and right-align the rest (numbers & percents).
            // This is just a style choice; adjust as desired.
            for (size_t c = 0; c < row.size(); ++c) {
                if (c == 0) {
                    // Left align the first column
                    std::cout << std::left << std::setw(colWidths[c]) << row[c];
                } else {
                    // Right align numeric columns
                    std::cout << std::right << std::setw(colWidths[c]) << row[c];
                }
                // Add a space after each column (except maybe the last)
                if (c + 1 < row.size()) {
                    std::cout << "  ";
                }
            }
            std::cout << "\n";
        }

        // Utility to print a separator line that matches the computed column widths
        void printSeparator(const std::vector<size_t> &colWidths) {
            // We'll match exactly the total text width, plus one space per gap between columns
            // (the same spacing used in printRow).
            // If you prefer a single continuous line, you can sum up everything and use '-' for the entire length.
            for (size_t c = 0; c < colWidths.size(); ++c) {
                // We add 2 for the two spaces after columns (except for the last column).
                // But let's do a simpler approach: just use '-' repeated colWidth times,
                // plus 2 for the spacing. For the last column we might skip the extra 2.
                std::cout << std::string(colWidths[c], '-');
                if (c + 1 < colWidths.size()) {
                    std::cout << "  ";
                }
            }
            std::cout << "\n";
        }
    };
}

#endif //HEIPROMAP_SMALL_STATISTIC_COLLECTOR_H
