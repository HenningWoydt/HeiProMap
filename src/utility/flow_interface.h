#ifndef HEIPROMAP_FLOW_INTERFACE_H
#define HEIPROMAP_FLOW_INTERFACE_H

#include <vector>
#include <memory>
#include "../definitions.h"
#include "macros.h"
#include "flow.h"

namespace HeiProMap {

template <typename captype, typename tcaptype, typename flowtype>
class IFlowAlgorithm {
public:
    virtual ~IFlowAlgorithm() = default;

    virtual void initialize(size_t n) = 0;
    virtual void add(vertex_t u, vertex_t v, weight_t w) = 0;
    virtual void add_s_edge(vertex_t v, weight_t w) = 0;
    virtual void add_t_edge(vertex_t v, weight_t w) = 0;
    virtual void solve() = 0;
    virtual void get_cut(std::vector<u8> &is_left) = 0;
    virtual void build_residual_network(ResidualFlowNetwork &residual_g) = 0;
    virtual void print() const = 0;
};

} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_INTERFACE_H
