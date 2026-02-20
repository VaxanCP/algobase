#pragma once

#include <cassert>
#include <vector>

#include "../internal/base/typing.hpp"
#include "./graph.hpp"
/*
@class/graph.hpp
@internal/base/typing.hpp
*/
// makecode
namespace algo {
class heavy_lignt_tree {
  struct hld_node {
    int parent;
    int depth;
    int head;
    int time;
  };
  using node_type = hld_node;
public:
  heavy_lignt_tree() noexcept = default;
  template <typename Weight>
  explicit heavy_lignt_tree(const adjacency_list<Weight> &adj,
                            int root = 0) noexcept
      : nodes_(adj.num_nodes()),
        size_(adj.num_nodes()),
        order_(adj.num_nodes()) {
#if !defined(NDEBUG)
    assert(0 <= root && root < num_nodes());
#endif
    decompose(adj, root);
  }
  int lca(int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        u = parent(head(u));
      } else {
        v = parent(head(v));
      }
    }
    return depth(u) < depth(v) ? u : v;
  }
  int dist(int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    return depth(u) + depth(v) - 2 * depth(lca(u, v));
  }
  int head(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return nodes_[u].head;
  }
  int time(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return nodes_[u].time;
  }
  int leave(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return time(u) + size(u);
  }
  int depth(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return nodes_[u].depth;
  }
  int parent(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return nodes_[u].parent;
  }
  int size(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return size_[u];
  }
  int node_index(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return time(u);
  }
  int edge_index(int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(v == parent(u) || u == parent(v));
#endif
    return parent(u) == v ? time(u) : time(v);
  }
  int order(int i) const {
#if !defined(NDEBUG)
    assert(0 <= i && i < num_nodes());
#endif
    return order_[i];
  }
  /**
   * \brief Return kth ancestor of u
   *
   * \param u
   * \param k
   * \return
   */
  int ancestor(int u, int k) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(k >= 0);
#endif
    if (depth(u) < k) { return null(); }
    while (depth(u) - depth(head(u)) < k) {
      k -= depth(u) - depth(head(u)) + 1;
      u = parent(head(u));
    }
    return order(time(u) - k);
  }
  /**
   * \brief Return the vertex on the path from u to v which is k distance
   * away from u
   *
   * \param u
   * \param v
   * \param k
   * \return
   */
  int move_to(int u, int v, int k) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(k >= 0);
#endif
    const int w = lca(u, v);
    const int uw = depth(u) - depth(w);
    const int vw = depth(v) - depth(w);
    if (k <= uw) { return ancestor(u, k); }
    return k <= uw + vw ? ancestor(v, uw + vw - k) : null();
  }
  /**
   * \brief Check if v is ancestor of u
   *
   * \param u
   * \param v
   * \return
   * \return
   */
  bool is_ancestor_of(int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    return time(v) <= time(u) && time(u) < leave(v);
  }
  /**
   * \brief Check if x is on the path between u and v
   *
   * \param x
   * \param u
   * \param v
   * \return
   * \return
   */
  bool on_path(int x, int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(0 <= x && x < num_nodes());
#endif
    return is_ancestor_of(x, lca(u, v)) &&
           (is_ancestor_of(u, x) || is_ancestor_of(v, x));
  }
  /**
   * \brief Get path between u and v
   *
   * \param u
   * \param v
   * \return
   */
  std::vector<int> get_path(int u, int v) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    const int w = lca(u, v);
    const int uw = depth(u) - depth(w);
    const int vw = depth(v) - depth(w);
    const int uv = uw + vw;
    std::vector<int> path(uv + 1);
    for (int i = 0; u != w; u = parent(u), ++i) { path[i] = u; }
    path[uw] = w;
    for (int i = uv; v != w; v = parent(v), --i) { path[i] = v; }
    return path;
  }
  /**
   * \brief Return vertex segments.
   *
   * \param u
   * \param v
   * \return
   */
  auto vertex_segments(int u, int v) const -> std::vector<std::array<int, 2>> {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    std::vector<std::array<int, 2>> segments{};
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        segments.push_back({time(head(u)), time(u) + 1});
        u = parent(head(u));
      } else {
        segments.push_back({time(head(v)), time(v) + 1});
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      segments.push_back({time(v), time(u) + 1});
    } else {
      segments.push_back({time(u), time(v) + 1});
    }
    return segments;
  }
  /**
   * \brief Return oriented vertex segments. Let w = lca(u,v), then for u to w,
   * r to l and w to v, l to r
   *
   * \param u
   * \param v
   * \return
   */
  auto oriented_vertex_segments(int u, int v) const
      -> std::array<std::vector<std::array<int, 2>>, 2> {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    std::vector<std::array<int, 2>> up{}, down{};
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        up.push_back({time(head(u)), time(u) + 1});
        u = parent(head(u));
      } else {
        down.push_back({time(head(v)), time(v) + 1});
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      up.push_back({time(v), time(u) + 1});
    } else {
      down.push_back({time(u), time(v) + 1});
    }
    std::reverse(down.begin(), down.end());
    return {std::move(up), std::move(down)};
  }
  /**
   * \brief Visit vertex segments. Let w = lca(u,v), and suppose that the
   * vertices which lies on the path between u to v will lies on the union of
   * disjoint segments of the form [L,R). Then fn will be called with all these
   * segments. Note that order of which segments will be called is undetermined.
   *
   * \tparam Fn
   * \param u
   * \param v
   * \param fn
   */
  template <std::invocable<int, int> Fn>
  void visit_vertex_segments(int u, int v, const Fn &fn) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        fn(time(head(u)), time(u) + 1);
        u = parent(head(u));
      } else {
        fn(time(head(v)), time(v) + 1);
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      fn(time(v), time(u) + 1);
    } else {
      fn(time(u), time(v) + 1);
    }
  }
  /**
   * \brief Visit vertex segments. Let w = lca(u,v). For the path from u to w,
   * suppose the path is decomposed into the union of disjoint segments in the
   * form [L0, R0), [L1, R1)... [Lk, Rk) such that going through from Ri-1
   * to Li for each segments in the forward order gives the sequence
   * of vertices which lies on the path between u to w. Then f1 will be called
   * with [L0, R0), [L1, R1),... [Lk, Rk) in this order. Similary, for the path
   * from w to v, suppose tha path is decomposed into the union of disjoint
   * segments in the form [L0, R0), [L1, R1),... [Lk, Rk) such that going
   * through from Li to Ri-1 in reverse order of segments gives the sequence of
   * vertices which lies on the path from w to v. Then f2 will be called with
   * [L0,R0), [L1, R1)... [Lk, Rk) in this order.
   *
   * \tparam F1
   * \tparam F2
   * \param u
   * \param v
   * \param fn
   */
  template <std::invocable<int, int> F1, std::invocable<int, int> F2>
  void visit_vertex_segments(int u, int v, const F1 &f1, const F2 &f2) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        f1(time(head(u)), time(u) + 1);
        u = parent(head(u));
      } else {
        f2(time(head(v)), time(v) + 1);
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      f1(time(v), time(u) + 1);
    } else {
      f2(time(u), time(v) + 1);
    }
  }
  auto edge_segments(int u, int v) const -> std::vector<std::array<int, 2>> {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    std::vector<std::array<int, 2>> segments{};
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        segments.push_back({time(head(u)), time(u) + 1});
        u = parent(head(u));
      } else {
        segments.push_back({time(head(v)), time(v) + 1});
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      segments.push_back({time(v) + 1, time(u) + 1});
    } else if (depth(u) < depth(v)) {
      segments.push_back({time(u) + 1, time(v) + 1});
    }
    return segments;
  }
  auto oriented_edge_segments(int u, int v) const
      -> std::array<std::vector<std::array<int, 2>>, 2> {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    std::vector<std::array<int, 2>> up{}, down{};
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        up.push_back({time(head(u)), time(u) + 1});
        u = parent(head(u));
      } else {
        down.push_back({time(head(v)), time(v) + 1});
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      up.push_back({time(v) + 1, time(u) + 1});
    } else if (depth(u) < depth(v)) {
      down.push_back({time(u) + 1, time(v) + 1});
    }
    std::reverse(down.begin(), down.end());
    return {std::move(up), std::move(down)};
  }
  template <detail::function<void(int, int)> Fn>
  void visit_edge_segments(int u, int v, const Fn &fn) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        fn(time(head(u)), time(u) + 1);
        u = parent(head(u));
      } else {
        fn(time(head(v)), time(v) + 1);
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      fn(time(v) + 1, time(u) + 1);
    } else if (depth(u) < depth(v)) {
      fn(time(u) + 1, time(v) + 1);
    }
  }
  template <detail::function<void(int, int)> F1,
            detail::function<void(int, int)> F2>
  void visit_edge_segments(int u, int v, const F1 &f1, const F2 &f2) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    while (head(u) != head(v)) {
      if (depth(head(u)) > depth(head(v))) {
        f1(time(head(u)), time(u) + 1);
        u = parent(head(u));
      } else {
        f2(time(head(v)), time(v) + 1);
        v = parent(head(v));
      }
    }
    if (depth(u) > depth(v)) {
      f1(time(v) + 1, time(u) + 1);
    } else if (depth(u) < depth(v)) {
      f2(time(u) + 1, time(v) + 1);
    }
  }
  std::vector<int> compressed() const { return order_; }
  int num_nodes() const { return static_cast<int>(order_.size()); }
  int null() const { return Null - 1; }
private:
  template <typename Weight>
  void decompose(const adjacency_list<Weight> &adj, int root = 0) {
    std::vector<int> heavy(num_nodes(), Null);
    std::fill(size_.begin(), size_.end(), 1);
    const auto dfs_tree = [&](auto self, int u, int p) -> void {
      int big_child = Null, max_subtree = 0;
      for (const int v : adj.neighbors(u)) {
        if (v == p) { continue; }
        nodes_[v].parent = u;
        nodes_[v].depth = nodes_[u].depth + 1;
        self(self, v, u);
        size_[u] += size_[v];
        if (max_subtree < size_[v]) {
          max_subtree = size_[v];
          big_child = v;
        }
      }
      heavy[u] = big_child;
    };
    int timer = 0;
    const auto dfs_hld = [&](auto self, int u, int h) -> void {
      nodes_[u].head = h;
      nodes_[u].time = timer;
      order_[timer++] = u;
      if (heavy[u] != Null) { self(self, heavy[u], h); }
      for (const int v : adj.neighbors(u)) {
        if (v != nodes_[u].parent && v != heavy[u]) { self(self, v, v); }
      }
    };
    nodes_[root].parent = Null;
    dfs_tree(dfs_tree, root, Null);
    dfs_hld(dfs_hld, root, root);
#if !defined(NDEBUG)
    assert(timer == num_nodes());
#endif
  }
  static constexpr int Null = -1;
  std::vector<node_type> nodes_;
  std::vector<int> size_;
  std::vector<int> order_;
};
}  // namespace algo