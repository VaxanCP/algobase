#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <span>
#include <vector>

#include "../internal/base/typing.hpp"
#include "./graph.hpp"
#include "./queue.hpp"

/*
@internal/base/typing.hpp
@class/queue.hpp
@class/graph.hpp
*/

// makecode
namespace algo {
template <detail::arithmetic Flow>
class maximum_flow {
  struct flow_edge_map {
    int v;
    int rev;
    Flow cap;
  };
  struct flow_node_map {
    int head;
    int level;
  };
  static constexpr int Undefined = -1;
  static constexpr Flow FlowInf = std::numeric_limits<Flow>::max() / 2;
public:
  using edge_type = edge<Flow>;
  maximum_flow() noexcept = default;
  explicit maximum_flow(int n) noexcept : num_nodes_{n}, edges_{} {}
  int add_edge(int u, int v, Flow cap) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(0 <= cap);
#endif
    const int id = num_arcs();
    edges_.emplace_back(u, v, cap);
    return id;
  }
  void set_capacity(int e, Flow new_cap) {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
    assert(0 <= new_cap);
#endif
    edges_[e].w = new_cap;
  }
  void add_capacity(int e, Flow d) {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    set_capacity(e, get_capacity(e) + d);
  }
  Flow get_capacity(int e) const {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    return edges_[e].w;
  }
  const edge_type& get_edge(int e) const {
#if !defined(NDEBUG)
    assert(0 <= e && e < num_arcs());
#endif
    return edges_[e];
  }
  /**
   * \brief Compute the maximum flow of the network.
   *
   * \param s
   * \param t
   * \return
   */
  Flow max_flow(int s, int t) const {
#if !defined(NDEBUG)
    assert(0 <= s && s < num_nodes());
    assert(0 <= t && t < num_nodes());
    assert(s != t);
#endif
    // This includes reversed edge.
    std::vector<flow_edge_map> elist(2 * num_arcs());
    std::vector<int> start(num_nodes() + 1);
    std::vector<int> edge_id(num_arcs());
    std::vector<flow_node_map> nodes(num_nodes());
    // This includes reversed edge
    std::vector<int> next_edge(2 * num_arcs());
    prepare(elist, start, edge_id);
    Flow total_flow = 0;
    while (reachable(s, t, start, elist, next_edge, nodes)) {
      total_flow += push_flow(s, t, FlowInf, next_edge, elist, nodes);
    }
    return total_flow;
  }
  /**
   * \brief Compute the maximum flow of the network and store one of assignment
   * of flows into to provided elements in range [flow.begin(), flow.end())
   *
   * \param s
   * \param t
   * \param flow
   * \return Flow
   */
  Flow max_flow(int s, int t, std::span<Flow> flow) const {
#if !defined(NDEBUG)
    assert(0 <= s && s < num_nodes());
    assert(0 <= t && t < num_nodes());
    assert(s != t);
    assert(flow.size() <= edges_.size());
#endif
    // This includes reversed edge.
    std::vector<flow_edge_map> elist(2 * num_arcs());
    std::vector<int> start(num_nodes() + 1);
    std::vector<int> edge_id(num_arcs());
    std::vector<flow_node_map> nodes(num_nodes());
    // This includes reversed edge
    std::vector<int> next_edge(2 * num_arcs());
    prepare(elist, start, edge_id);
    Flow total_flow = 0;
    while (reachable(s, t, start, elist, next_edge, nodes)) {
      total_flow += push_flow(s, t, FlowInf, next_edge, elist, nodes);
    }
    for (size_t i = 0; i < flow.size(); ++i) {
      flow[i] = edges_[i].w - elist[edge_id[i]].cap;
    }
    return total_flow;
  }
  /**
   * \brief Compute the minimum s-t cut of the network and return one of the
   * assignments.
   *
   * \param s
   * \param t
   * \return One of possible assignments. ith elements indicate that whether the
   * vertex i is in T side or not. (false means the vertex is in S side,
   * otherwise true).
   */
  std::vector<bool> min_cut(int s, int t) const {
#if !defined(NDEBUG)
    assert(0 <= s && s < num_nodes());
    assert(0 <= t && t < num_nodes());
    assert(s != t);
#endif
    // This includes reversed edge.
    std::vector<flow_edge_map> elist(2 * num_arcs());
    std::vector<int> start(num_nodes() + 1);
    std::vector<int> edge_id(num_arcs());
    std::vector<flow_node_map> nodes(num_nodes());
    // This includes reversed edge
    std::vector<int> next_edge(2 * num_arcs());
    prepare(elist, start, edge_id);
    Flow total_flow = 0;
    while (reachable(s, t, start, elist, next_edge, nodes)) {
      total_flow += push_flow(s, t, FlowInf, next_edge, elist, nodes);
    }
    std::vector<bool> right(num_nodes());
    std::transform(nodes.begin(), nodes.end(), right.begin(),
                   [](flow_node_map nd) { return nd.level == Undefined; });
    return right;
  }
  /**
   * \brief Return the minimum s-t cut.
   *
   * \param s
   * \param t
   * \return
   */
  std::vector<bool> min_cut(int s, int t, std::span<Flow> flow) const {
#if !defined(NDEBUG)
    assert(0 <= s && s < num_nodes());
    assert(0 <= t && t < num_nodes());
    assert(s != t);
    assert(flow.size() <= edges_.size());
#endif
    // This includes reversed edge.
    std::vector<flow_edge_map> elist(2 * num_arcs());
    std::vector<int> start(num_nodes() + 1);
    std::vector<int> edge_id(num_arcs());
    std::vector<flow_node_map> nodes(num_nodes());
    // This includes reversed edge
    std::vector<int> next_edge(2 * num_arcs());
    prepare(elist, start, edge_id);
    Flow total_flow = 0;
    while (reachable(s, t, start, elist, next_edge, nodes)) {
      total_flow += push_flow(s, t, FlowInf, next_edge, elist, nodes);
    }
    for (size_t i = 0; i < flow.size(); ++i) {
      flow[i] = edges_[i].w - elist[edge_id[i]].cap;
    }
    std::vector<bool> right(num_nodes());
    std::transform(nodes.begin(), nodes.end(), right.begin(),
                   [](flow_node_map nd) { return nd.level == Undefined; });
    return right;
  }
  int num_nodes() const { return num_nodes_; }
  int num_arcs() const { return static_cast<int>(edges_.size()); }
private:
  /**
   * \brief Construct residual network.
   *
   * \param elist
   * \param start
   * \param edge_id
   */
  void prepare(std::vector<flow_edge_map>& elist, std::vector<int>& start,
               std::vector<int>& edge_id) const {
    for (const auto& [u, v, cap] : edges_) {
      start[u]++;
      start[v]++;
    }
    for (int u = 1; u <= num_nodes(); ++u) { start[u] += start[u - 1]; }
    for (int i = num_arcs() - 1; i >= 0; --i) {
      const auto& [u, v, cap] = edges_[i];
      const int e1 = --start[u];
      const int e0 = --start[v];
      elist[e1] = {.v = v, .rev = e0, .cap = cap};
      elist[e0] = {.v = u, .rev = e1, .cap = 0};
      edge_id[i] = e1;
    }
  }
  bool reachable(int s, int t, const std::vector<int>& start,
                 const std::vector<flow_edge_map>& elist,
                 std::vector<int>& next_edge,
                 std::vector<flow_node_map>& nodes) const {
    std::fill(next_edge.begin(), next_edge.end(), Undefined);
    std::fill(nodes.begin(), nodes.end(), flow_node_map{Undefined, Undefined});
    buffer_queue<int> que;
    que.push(s);
    nodes[s].level = 0;
    while (!que.empty()) {
      const int u = que.front();
      que.pop();
      for (int e = start[u]; e < start[u + 1]; ++e) {
        if (elist[e].cap == 0) { continue; }
        const int v = elist[e].v;
        if (nodes[v].level == Undefined) {
          nodes[v].level = nodes[u].level + 1;
          que.push(v);
        }
        if (nodes[u].level < nodes[v].level) {
          next_edge[e] = nodes[u].head;
          nodes[u].head = e;
        }
      }
    }
    return nodes[t].level != Undefined;
  }
  Flow push_flow(int u, int t, Flow flow_to_push,
                 const std::vector<int>& next_edge,
                 std::vector<flow_edge_map>& elist,
                 std::vector<flow_node_map>& nodes) const {
    if (u == t) { return flow_to_push; }
    Flow pushed_flow = 0;
    int& e = nodes[u].head;
    while (e != Undefined && flow_to_push > 0) {
      const int v = elist[e].v;
      const Flow next_flow = std::min(flow_to_push, elist[e].cap);
      const Flow improved = push_flow(v, t, next_flow, next_edge, elist, nodes);
      if (improved != 0) {
        elist[e].cap -= improved;
        elist[elist[e].rev].cap += improved;
        pushed_flow += improved;
        flow_to_push -= improved;
      }
      if (elist[e].cap == 0 || flow_to_push > 0) { e = next_edge[e]; }
    }
    return pushed_flow;
  }
  int num_nodes_;
  std::vector<edge_type> edges_;
};
}  // namespace algo
