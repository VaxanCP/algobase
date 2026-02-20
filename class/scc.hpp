#pragma once

#include <cassert>
#include <vector>

#include "../internal/base/def.hpp"
#include "./graph.hpp"
/*
@class/graph.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo {
/**
 * \brief Strongly connected components algorithm implementation.
 * Precomputation is done in constructor.
 *
 */
class scc_graph {
private:
  static constexpr int Undefined = -1;
  struct scc_node {
    int time = Undefined;
    int low = Undefined;
    int comp = Undefined;
  };
  using node_type = scc_node;
public:
  using graph_type = adjacency_list<>;
  scc_graph() noexcept = default;
  explicit scc_graph(const graph_type& adj) noexcept
      : num_nodes_{0}, num_cc_{0}, nodes_(adj.num_nodes()) {
    build_scc(adj);
  }
  int id(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return nodes_[u].comp;
  }
  int num_nodes() const { return num_nodes_; }
  int num_components() const { return num_cc_; }
private:
  /**
   * \brief This is only called on construction
   *
   * \param adj graph
   */
  void build_scc(const graph_type& adj) {
    std::vector<int> rec_stack;
    rec_stack.reserve(adj.num_nodes());
    const auto dfs = [&](auto dfs, int u) -> void {
      nodes_[u].time = nodes_[u].low = num_nodes_++;
      rec_stack.emplace_back(u);
      for (const int v : adj.neighbors(u)) {
        if (nodes_[v].time == Undefined) {
          // v is not visited
          dfs(dfs, v);
          nodes_[u].low = std::min(nodes_[u].low, nodes_[v].low);
        } else if (nodes_[v].comp == Undefined) {
          // v is on stack
          nodes_[u].low = std::min(nodes_[u].low, nodes_[v].time);
        }
      }
      if (nodes_[u].low == nodes_[u].time) {
        // Found a scc
        while (rec_stack.back() != u) {
          const int w = rec_stack.back();
          nodes_[w].comp = num_cc_;
          rec_stack.pop_back();
        }
        rec_stack.pop_back();
        nodes_[u].comp = num_cc_++;
      }
    };
    for (int u = 0; u < adj.num_nodes(); ++u) {
      if (nodes_[u].time == Undefined) { dfs(dfs, u); }
    }
  }
  int num_nodes_;
  int num_cc_;
  std::vector<node_type> nodes_;
};
}  // namespace algo