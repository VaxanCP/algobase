#pragma once

#include <cassert>
#include <optional>
#include <vector>

#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
#include "./maxflow.hpp"

/*
@class/maxflow.hpp
@internal/base/typing.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo {
template <detail::arithmetic Flow>
class circulation {
  struct flow_edge {
    int u;
    int v;
    Flow low;
    Flow cap;
  };
public:
  using edge_type = flow_edge;
  using solver_type = maximum_flow<Flow>;
  circulation() noexcept = default;
  explicit circulation(int n) noexcept : supply_(n), edges_{} {}
  int add_edge(int u, int v, Flow low, Flow cap) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(0 <= low && low <= cap);
#endif
    const int id = num_arcs();
    edges_.emplace_back(u, v, low, cap);
    return id;
  }
  void set_lower(int e, Flow new_low) {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
    assert(0 <= new_low && get_capacity(e));
#endif
    edges_[e].low = new_low;
  }
  void set_capacity(int e, Flow new_cap) {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
    assert(0 <= new_cap && get_lower(e) <= new_cap);
#endif
    edges_[e].cap = new_cap;
  }
  void set_supply(int u, Flow new_supply) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    supply_[u] = new_supply;
  }
  void set_demand(int u, Flow new_demand) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    supply_[u] = -new_demand;
  }
  void add_demand(int u, Flow d) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    supply_[u] -= d;
  }
  void add_supply(int u, Flow d) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    supply_[u] += d;
  }
  Flow get_lower(int e) const {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    return edges_[e].low;
  }
  Flow get_capacity(int e) const {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    return edges_[e].cap;
  }
  Flow get_supply(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return supply_[u];
  }
  Flow get_demand(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return -get_supply(u);
  }
  const flow_edge& get_edge(int e) const {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    return edges_[e];
  }
  /**
   * \brief Check if the circulation network is feasible or not.
   *
   * \return true if it's possible to assign the flows to each edges to satisfy
   * the demands requirements for every vertex.
   * \return false Otherwise
   */
  bool is_feasible() const {
    const int s_dummy = num_nodes();
    const int t_dummy = num_nodes() + 1;
    solver_type mf_net(t_dummy + 1);
    Flow total_supply = 0;
    Flow total_demand = 0;
    std::vector<Flow> excess = supply_;
    for (const auto& [u, v, low, cap] : edges_) {
      excess[u] -= low;
      excess[v] += low;
      mf_net.add_edge(u, v, cap - low);
    }
    for (int u = 0; u < num_nodes(); ++u) {
      if (excess[u] > 0) {
        mf_net.add_edge(s_dummy, u, excess[u]);
        total_supply += excess[u];
      } else if (excess[u] < 0) {
        mf_net.add_edge(u, t_dummy, -excess[u]);
        total_demand += -excess[u];
      }
    }
    const Flow total_flow = mf_net.max_flow(s_dummy, t_dummy);
    return total_flow == total_demand && total_flow == total_supply;
  }
  /**
   * \brief Check if the circulation network is feasible and if so, returns one
   * of possible assignments of flows.
   *
   * \return One of possible assignments of flows or nullopt if it's not
   * feasible.
   */
  auto feasible_circulation() const -> std::optional<std::vector<Flow>> {
    const int s_dummy = num_nodes();
    const int t_dummy = num_nodes() + 1;
    solver_type mf_net(t_dummy + 1);
    Flow total_supply = 0;
    Flow total_demand = 0;
    std::vector<Flow> excess = supply_;
    for (const auto& [u, v, low, cap] : edges_) {
      excess[u] -= low;
      excess[v] += low;
      mf_net.add_edge(u, v, cap - low);
    }
    for (int u = 0; u < num_nodes(); ++u) {
      if (excess[u] > 0) {
        mf_net.add_edge(s_dummy, u, excess[u]);
        total_supply += excess[u];
      } else if (excess[u] < 0) {
        mf_net.add_edge(u, t_dummy, -excess[u]);
        total_demand += -excess[u];
      }
    }
    // Compute assignment
    std::vector<Flow> flows(num_arcs());
    const Flow total_flow = mf_net.max_flow(s_dummy, t_dummy, flows);
    if (total_flow != total_demand || total_flow != total_supply) {
      return std::nullopt;
    }
    return {std::move(flows)};
  }
  int num_nodes() const { return static_cast<int>(supply_.size()); }
  int num_arcs() const { return static_cast<int>(edges_.size()); }
private:
  std::vector<Flow> supply_;
  std::vector<edge_type> edges_;
};
}  // namespace algo