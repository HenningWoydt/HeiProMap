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

#ifndef HEIPROMAP_DEEP_REBALANCER_H
#define HEIPROMAP_DEEP_REBALANCER_H

#include "../../serial/serial_definitions_1.h"
#include "../../serial/serial_definitions_2.h"
#include "../../serial/serial_definitions_3.h"
#include "../../serial/utility/qap.h"

namespace HeiProMap {
    /**
     * As described in 4.1
     * https://pdf.sciencedirectassets.com/272438/1-s2.0-S0743731500X01062/1-s2.0-S0743731597914106/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIAZnp1PThtso0fo2rZlDvbPcyqdj2PJ910oo28Dabtu8AiEApnrHuGgdRHHp%2Ffqy%2BibX%2BBmdiGXm0n7bc7rAclZJO9sqvAUIkf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDKDDKA48BcHpa4I5jSqQBYkyCqhznzJZkNwKJ1wwpnsEUYPGR5so4QtPjkDWFvE5bF9cxqnLMczXiC%2F8fCLOL4WQMpNrL55K1C%2Fhwm6SjvegDr1byeuabeu5sjtiXW4Di%2BBpBVxPWODBhz2HjZndDdvdPv1RL%2BYiRONFN8gIcAPZFXSXSzKJJ5WkMCrNA9C2YVgex2hEwsE5fbIkUVvYHhtzEtUALJW2Aq7XyP%2Fg5jEBFeDMFFleqgPpATQ78ATwWQ6htaFt%2F0UBEQC1OsvfbOsJPpUJo%2BV2iUUflKQSH84X4ho%2FKGXKPtrAlum6qbbweBMjK9j33SHgZUGIYUemKk28sb2GuFh98zh0D4Fb94qr8Py2n5eo%2FLaePUQEY3mBIBL7FnaEm%2F%2FZEcJ1aO5OS4juyRWo7IIWro0w4jYMymvycQMsFThl7BT4XtRZWGX4MaUHOVW8k7cnBPA21vJpUc6xA%2FYxxbZFQMbQxq1CiHz0UYM50fSp3HZvs%2FJ0Ju4V1%2FkSnNE7ZhMHcfIoODKoZXbjIUBoxRdgCKYJiopNEnNwJKuhqQD4cC0G3LGJ9cRwHg2OvIerNFpaIZ8OQHh1hkxQwZx2zAxWSLa%2F0Ac%2BD44n8T7h%2F3zpJhZIXnPyTEdR%2B%2BZ%2FnuwJYcTYmPpHb89rsrLCqcYBkWs%2Fc9no4%2F6cAaqZwxV3sfi8S%2FQFQrxXGoI4fLLRZyeElCJnxYZy2RTZhzgagtLJRC8nmXmzobu7ouOD%2FAn8NgrV9LmzL199IrlnYyj%2FvrYBz8QL8lhbEj8wClBfsJ%2F6Gx07ASXQtwTfgYXLVLQxap3lgebev4wcajFhHrSCp9VnbWUbevUMTrBQDRhy3kHMJJKC1tcK%2BhDQCfHghSdTnJqLW9hEsTk1ZHgRMN27y8IGOrEB3bJJ%2FfqeyGvYkPdDucmUxeeBmWjp84aHXmcswSwvGzKXyIf%2FxxjBtIcQgfoKaq8Vnm%2Fwh08bleRv%2BDoIKkZSN6X5xztg8EAioglSBaPMjnWDDXRr6FTTOLrHpAmGC0iYuJwPURGG6fZcxyr%2F0ElH3pMaEeZu3JfiekM0YzNksRa6ps8Keuq51kO%2FUS5QQiUER46q52KaIxP7ftQfS3Uyi77CdauKHl7YCs1QlAl3CZWT&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250618T161716Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY54YD5F37%2F20250618%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e0857bbffd3e3c7ec5dc454b245b5d11d8e212fca2fbe064d4d59378d5c74015&hash=1230f7c24327de6d67982e0966b74873ac118f60185fe7f3c2c85cd62b1e2d76&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0743731597914106&tid=spdf-c5bba773-ff59-45b7-b4df-60a177db705e&sid=0d1e732f9ff0a745948961e7b34832015f6dgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1e035d51570c555b560a52&rr=951c188d7c8e2c0c&cc=de
     */
    class DeepRebalancer {
        RandomEngine* random_engine = nullptr;

    public:
        void initialize(RandomEngine& t_random_engine) {
            random_engine = &t_random_engine;
        }

        void rebalance(const deep_graph_t& g,
                       deep_p_manager_t& p_manager,
                       deep_bv_manager_t& bv_manager,
                       deep_q_graph_t& q_graph,
                       deep_d_oracle_t& d_oracle,
                       partition_t k) {
            std::vector<vertex_t> boundary;

            while (true) {
                // collect all vertices
                boundary.clear();
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                        forall_bv_id_iu(bv_manager, id, i, u)
                            {
                                boundary.push_back(u);
                            }
                        endfor
                    }
                }

                // shuffle them
                std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                bool made_move = false;
                for (vertex_t u : boundary) {
                    partition_t u_id  = p_manager[u];
                    weight_t u_weight = g.weight(u);

                    if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }

                    // move this vertex such that the other block is not overloaded, and we have the smallest loss
                    partition_t best_id = 0;
                    s64 best_qap_delta  = -std::numeric_limits<s64>::max();
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];

                            if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_id        = v_id;
                            }
                        }
                    endfor

                    for (partition_t v_id : q_graph.lowest_level_neighborhood(u_id)) {
                        if (!p_manager.is_active(v_id)) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta) {
                            best_qap_delta = qap_delta;
                            best_id        = v_id;
                        }
                    }

                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        p_manager.move(u, u_weight, u_id, best_id);

                        made_move = true;
                    }
                }
                if (!made_move) {
                    break;
                }
            }
        }

        void rebalance_last_layer(const deep_graph_t& g,
                                  deep_p_manager_t& p_manager,
                                  deep_bv_manager_t& bv_manager,
                                  deep_q_graph_t& q_graph,
                                  deep_d_oracle_t& d_oracle,
                                  partition_t k) {
            rebalance(g, p_manager, bv_manager, q_graph, d_oracle, k);

            std::vector<vertex_t> boundary;

            bool global_move = true;
            while (global_move) {
                global_move = false;
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) <= p_manager.get_lmax(id)) { continue; }
                    // collect all vertices
                    boundary.clear();
                    forall_bv_id_iu(bv_manager, id, i, u)
                        {
                            boundary.push_back(u);
                        }
                    endfor

                    // shuffle them
                    std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                    for (vertex_t u : boundary) {
                        partition_t u_id  = p_manager[u];
                        if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }
                        partition_t best_id = 0;
                        s64 best_qap_delta  = -std::numeric_limits<s64>::max();
                        weight_t u_weight   = g.weight(u);
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];

                                if (v_id == id) { continue; }
                                if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);

                                if (qap_delta > best_qap_delta) {
                                    best_qap_delta = qap_delta;
                                    best_id        = v_id;
                                }
                            }
                        endfor

                        for (partition_t v_id : q_graph.lowest_level_neighborhood(u_id)) {
                            if (v_id == id) { continue; }
                            if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);

                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_id        = v_id;
                            }
                        }

                        if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);

                            global_move = true;
                        }
                    }
                }
            }

            // No global move made. If still overloaded, make desperate moves
            bool overloaded = true;
            while (overloaded) {
                overloaded = false;

                // collect all vertices
                boundary.clear();
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                        overloaded = true;
                        forall_bv_id_iu(bv_manager, id, i, u)
                            {
                                boundary.push_back(u);
                            }
                        endfor
                    }
                }

                // shuffle them
                std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                for (vertex_t u : boundary) {
                    partition_t u_id = p_manager[u];
                    if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }

                    partition_t best_id = 0;
                    s64 best_qap_delta  = -std::numeric_limits<s64>::max();
                    weight_t u_weight   = g.weight(u);

                    for (partition_t v_id = 0; v_id < k; ++v_id) {
                        if (v_id == u_id) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                        if (qap_delta > best_qap_delta) {
                            best_qap_delta = qap_delta;
                            best_id        = v_id;
                        }
                    }

                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        p_manager.move(u, u_weight, u_id, best_id);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_REBALANCER_H
