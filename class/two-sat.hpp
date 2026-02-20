#pragma once

#include <cassert>
#include <utility>
#include <vector>

#include "../internal/base/def.hpp"
#include "./graph.hpp"
#include "./scc.hpp"

/*
@internal/base/def.hpp
@class/scc.hpp
@class/graph.hpp
*/
// makecode
namespace algo {
class two_sat {
  using graph_type = scc_graph::graph_type;
  using edge_type = graph_type::edge_type;
public:
  two_sat() noexcept = default;
  explicit two_sat(int n) noexcept : num_nodes_{n}, edges_{} {}
  /**
   * \brief Add disjunction clause such that u OR v = true
   *
   * \param u
   * \param v
   * \param nu
   * \param nv
   */
  void add_disjunction(int u, int v, bool nu, bool nv) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    // If nu is true, 2u+1
    u = (2 * u) ^ nu;
    // If nv is true
    v = (2 * v) ^ nv;
    edges_.emplace_back(u ^ 1, v);
    edges_.emplace_back(v ^ 1, u);
  }
  /**
   * \brief Check if given CNF can be satisfiable
   *
   * \return
   * \return
   */
  bool is_feasible() const {
    // For each u, the node 2u represents the case when the term u is true, and
    // 2u+1 represents the case when the term u is false.
    // A CNF is feasible if and only if for each term u, 2u and 2u+1 do not
    // belong to the same strongly connected component of the implication graph.
    // The detailed explanation will be found here:
    // https://cp-algorithms.com/graph/2SAT.html
    const scc_graph scc(graph_type(directed, 2 * num_nodes(), edges_));
    for (int u = 0; u < num_nodes(); ++u) {
      if (scc.id(2 * u) == scc.id(2 * u + 1)) { return false; }
    }
    return true;
  }
  /**
   * \brief Solve the 2-sat and return if it is satisfiable or not.
   *
   * \return
   * \return
   */
  auto feasible_solution() const -> std::optional<std::vector<bool>> {
    // For each u, the node 2u represents the case when the term u is true, and
    // 2u+1 represents the case when the term u is false.
    // A CNF is feasible if and only if for each term u, 2u and 2u+1 do not
    // belong to the same strongly connected component of the implication graph.
    // For the construction, refer to: https://cp-algorithms.com/graph/2SAT.html
    const scc_graph scc(graph_type(directed, 2 * num_nodes(), edges_));
    std::vector<bool> ans(num_nodes());
    for (int u = 0; u < num_nodes(); ++u) {
      if (scc.id(2 * u) == scc.id(2 * u + 1)) { return std::nullopt; }
      ans[u] = scc.id(2 * u) < scc.id(2 * u + 1);
    }
    return {std::move(ans)};
  }
  int num_nodes() const { return num_nodes_; }
private:
  int num_nodes_;
  std::vector<edge_type> edges_;
};
}  // namespace algo