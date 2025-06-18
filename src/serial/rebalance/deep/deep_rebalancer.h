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

#include "../../serial_definitions_1.h"
#include "../../serial_definitions_2.h"
#include "../../serial_definitions_3.h"

namespace HeiProMap {
    /**
     * As described in
     * https://pdf.sciencedirectassets.com/272438/1-s2.0-S0743731500X01062/1-s2.0-S0743731597914106/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIAZnp1PThtso0fo2rZlDvbPcyqdj2PJ910oo28Dabtu8AiEApnrHuGgdRHHp%2Ffqy%2BibX%2BBmdiGXm0n7bc7rAclZJO9sqvAUIkf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDKDDKA48BcHpa4I5jSqQBYkyCqhznzJZkNwKJ1wwpnsEUYPGR5so4QtPjkDWFvE5bF9cxqnLMczXiC%2F8fCLOL4WQMpNrL55K1C%2Fhwm6SjvegDr1byeuabeu5sjtiXW4Di%2BBpBVxPWODBhz2HjZndDdvdPv1RL%2BYiRONFN8gIcAPZFXSXSzKJJ5WkMCrNA9C2YVgex2hEwsE5fbIkUVvYHhtzEtUALJW2Aq7XyP%2Fg5jEBFeDMFFleqgPpATQ78ATwWQ6htaFt%2F0UBEQC1OsvfbOsJPpUJo%2BV2iUUflKQSH84X4ho%2FKGXKPtrAlum6qbbweBMjK9j33SHgZUGIYUemKk28sb2GuFh98zh0D4Fb94qr8Py2n5eo%2FLaePUQEY3mBIBL7FnaEm%2F%2FZEcJ1aO5OS4juyRWo7IIWro0w4jYMymvycQMsFThl7BT4XtRZWGX4MaUHOVW8k7cnBPA21vJpUc6xA%2FYxxbZFQMbQxq1CiHz0UYM50fSp3HZvs%2FJ0Ju4V1%2FkSnNE7ZhMHcfIoODKoZXbjIUBoxRdgCKYJiopNEnNwJKuhqQD4cC0G3LGJ9cRwHg2OvIerNFpaIZ8OQHh1hkxQwZx2zAxWSLa%2F0Ac%2BD44n8T7h%2F3zpJhZIXnPyTEdR%2B%2BZ%2FnuwJYcTYmPpHb89rsrLCqcYBkWs%2Fc9no4%2F6cAaqZwxV3sfi8S%2FQFQrxXGoI4fLLRZyeElCJnxYZy2RTZhzgagtLJRC8nmXmzobu7ouOD%2FAn8NgrV9LmzL199IrlnYyj%2FvrYBz8QL8lhbEj8wClBfsJ%2F6Gx07ASXQtwTfgYXLVLQxap3lgebev4wcajFhHrSCp9VnbWUbevUMTrBQDRhy3kHMJJKC1tcK%2BhDQCfHghSdTnJqLW9hEsTk1ZHgRMN27y8IGOrEB3bJJ%2FfqeyGvYkPdDucmUxeeBmWjp84aHXmcswSwvGzKXyIf%2FxxjBtIcQgfoKaq8Vnm%2Fwh08bleRv%2BDoIKkZSN6X5xztg8EAioglSBaPMjnWDDXRr6FTTOLrHpAmGC0iYuJwPURGG6fZcxyr%2F0ElH3pMaEeZu3JfiekM0YzNksRa6ps8Keuq51kO%2FUS5QQiUER46q52KaIxP7ftQfS3Uyi77CdauKHl7YCs1QlAl3CZWT&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250618T161716Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY54YD5F37%2F20250618%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e0857bbffd3e3c7ec5dc454b245b5d11d8e212fca2fbe064d4d59378d5c74015&hash=1230f7c24327de6d67982e0966b74873ac118f60185fe7f3c2c85cd62b1e2d76&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0743731597914106&tid=spdf-c5bba773-ff59-45b7-b4df-60a177db705e&sid=0d1e732f9ff0a745948961e7b34832015f6dgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1e035d51570c555b560a52&rr=951c188d7c8e2c0c&cc=de
     */
    class DeepRebalancer {
        f64 time = 0;

    public:
        void rebalance(const graph_t &g,
                       deep_p_manager_t &p_manager,
                       deep_bv_manager_t &bv_manager,
                       deep_q_graph_t &q_graph,
                       partition_t k) {
            for (partition_t id = 0; id < k; ++id) {
                if (p_manager.is_active(id) && p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                    // move a boundary vertex that will underload the block, or at least reduce the weight

                    weight_t min_weight = p_manager.get_bweight(id) - p_manager.get_lmax(id);

                    vertex_t best_u = g.get_n();
                    weight_t best_weight = 0;
                    forall_bv_id_iu(bv_manager, id, i, u)
                        {
                            weight_t w = g.weight(u);

                            if (best_weight > min_weight) {
                                // look for a weight smaller than best weight, but larger equal min weight
                                if (w < best_weight && w >= min_weight) {
                                    best_u = u;
                                    best_weight = w;
                                }
                            } else {
                                if (w > best_weight) {
                                    best_u = u;
                                    best_weight = w;
                                }
                            }
                        }
                    endfor

                    if (best_u == g.get_n()) {
                        // no vertex found, continue
                        continue;
                    }

                    // first look in the neighborhood if we can move the vertex there
                    bool moved = false;
                    for (auto &[v_id, w]: q_graph.neighborhood(id)) {
                        if (p_manager.get_bweight(v_id) + best_weight <= p_manager.get_lmax(v_id)) {

                            bv_manager.move(g, p_manager, best_u, id, v_id);
                            q_graph.move(g, p_manager, best_u, id, v_id);
                            p_manager.move(best_u, best_weight, id, v_id);

                            moved = true;
                            break;
                        }
                    }

                    if (!moved) {
                        // if no place found checks all partitions
                        for (partition_t v_id = 0; v_id < k; ++v_id) {
                            if (p_manager.is_active(v_id) && p_manager.get_bweight(v_id) + best_weight <= p_manager.get_lmax(v_id)) {

                                bv_manager.move(g, p_manager, best_u, id, v_id);
                                q_graph.move(g, p_manager, best_u, id, v_id);
                                p_manager.move(best_u, best_weight, id, v_id);

                                moved = true;
                                break;
                            }
                        }
                    }

                    if (moved && p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                        id--; // check same partition again, if still overloaded
                    }

                    // the vertex cannot be moved, hopefully problem will solve itself after uncoarsening
                }
            }
        }


    };
}

#endif //HEIPROMAP_DEEP_REBALANCER_H
