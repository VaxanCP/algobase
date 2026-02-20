#pragma once

#include <algorithm>
#include <cassert>
#include <compare>
#include <ranges>
#include <span>
#include <vector>

#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
/*
@internal/base/typing.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo {
template <std::semiregular Weight = detail::monostate>
struct edge {
  int u;
  int v;
  Weight w;
  using arc_type = std::pair<int, Weight>;
};
template <>
struct edge<detail::monostate> {
  int u;
  int v;
  using arc_type = int;
};
/**
 * \brief Tag structs
 *
 */
struct directed_tag {};
struct undirected_tag {};

inline constexpr directed_tag directed{};
inline constexpr undirected_tag undirected{};
/**
 * \brief Adjacency list implementation
 *
 * \tparam Weight
 */
template <std::semiregular Weight = detail::monostate>
class adjacency_list {
public:
  using edge_type = edge<Weight>;
  using arc_type = typename edge_type::arc_type;
  adjacency_list() noexcept = default;
  adjacency_list(directed_tag, int n,
                 const std::vector<edge_type>& edges) noexcept;
  adjacency_list(undirected_tag, int n,
                 const std::vector<edge_type>& edges) noexcept;
  std::span<const arc_type> operator[](int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return {elist_.cbegin() + start_[u], elist_.cbegin() + start_[u + 1]};
  }
  std::span<arc_type> operator[](int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return {elist_.begin() + start_[u], elist_.begin() + start_[u + 1]};
  }
  auto neighbors(int u) const;
  int degree(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return start_[u + 1] - start_[u];
  }
  int num_nodes() const { return static_cast<int>(start_.size()) - 1; }
  int num_arcs() const { return static_cast<int>(elist_.size()); }
private:
  std::vector<int> start_;
  std::vector<arc_type> elist_;
};
template <std::semiregular Weight>
adjacency_list<Weight>::adjacency_list(
    directed_tag, int n, const std::vector<edge_type>& edges) noexcept
    : start_(n + 1), elist_(edges.size()) {
  for (const auto& e : edges) { start_[e.u]++; }
  for (int i = 1; i <= n; ++i) { start_[i] += start_[i - 1]; }
  for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
    elist_[--start_[it->u]] = {it->v, it->w};
  }
}
template <std::semiregular Weight>
adjacency_list<Weight>::adjacency_list(
    undirected_tag, int n, const std::vector<edge_type>& edges) noexcept
    : start_(n + 1), elist_(edges.size() * 2) {
  for (const auto& e : edges) {
    start_[e.u]++;
    start_[e.v]++;
  }
  for (int i = 1; i <= n; ++i) { start_[i] += start_[i - 1]; }
  for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
    elist_[--start_[it->u]] = {it->v, it->w};
    elist_[--start_[it->v]] = {it->u, it->w};
  }
}
template <>
adjacency_list<detail::monostate>::adjacency_list(
    directed_tag, int n, const std::vector<edge_type>& edges) noexcept
    : start_(n + 1), elist_(edges.size()) {
  for (const auto& e : edges) { start_[e.u]++; }
  for (int i = 1; i <= n; ++i) { start_[i] += start_[i - 1]; }
  for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
    elist_[--start_[it->u]] = {it->v};
  }
}
template <>
adjacency_list<detail::monostate>::adjacency_list(
    undirected_tag, int n, const std::vector<edge_type>& edges) noexcept
    : start_(n + 1), elist_(edges.size() * 2) {
  for (const auto& e : edges) {
    start_[e.u]++;
    start_[e.v]++;
  }
  for (int i = 1; i <= n; ++i) { start_[i] += start_[i - 1]; }
  for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
    elist_[--start_[it->u]] = {it->v};
    elist_[--start_[it->v]] = {it->u};
  }
}
template <std::semiregular Weight>
auto adjacency_list<Weight>::neighbors(int u) const {
#if !defined(NDEBUG)
  assert(0 <= u && u < num_nodes());
#endif
  return std::views::keys(operator[](u));
}
template <>
auto adjacency_list<detail::monostate>::neighbors(int u) const {
#if !defined(NDEBUG)
  assert(0 <= u && u < num_nodes());
#endif
  return operator[](u);
}
// Deduction Guides
template <typename Weight>
adjacency_list(directed_tag, int, const std::vector<edge<Weight>>&)
    -> adjacency_list<Weight>;
template <typename Weight>
adjacency_list(undirected_tag, int, const std::vector<edge<Weight>>&)
    -> adjacency_list<Weight>;
/**
 * \brief
 *
 * \tparam Weight
 * \param adj
 * \return std::vector<int>
 */
template <typename Weight>
std::vector<int> toposort(const adjacency_list<Weight>& adj) {
  std::vector<bool> vis(adj.num_nodes());
  std::vector<int> order;
  order.reserve(adj.num_nodes());
  const auto dfs = [&](auto dfs, int u) -> void {
    vis[u] = true;
    for (const int v : adj.neighbors(u)) {
      if (!vis[v]) { dfs(dfs, v); }
    }
    order.emplace_back(u);
  };
  for (int u = 0; u < adj.num_nodes(); ++u) {
    if (!vis[u]) { dfs(dfs, u); }
  }
  std::reverse(order.begin(), order.end());
  return order;
}
/**
 * \brief
 *
 * \tparam Weight
 * \param adj
 * \param root
 * \return std::vector<std::array<int, 2>>
 */
template <typename Weight>
auto preorder(const adjacency_list<Weight>& adj, int root = 0)
    -> std::vector<std::array<int, 2>> {
#if !defined(NDEBUG)
  assert(0 <= root && root < adj.num_nodes());
#endif
  std::vector<std::array<int, 2>> res(adj.num_nodes());
  int timer = 0;
  const auto dfs = [&](auto dfs, int u, int p) -> void {
    res[u][0] = timer++;
    for (const int v : adj.neighbors(u)) {
      if (v == p) { continue; }
      dfs(dfs, v, u);
    }
    res[u][1] = timer;
  };
  dfs(dfs, root, root);
  return res;
}
}  // namespace algo
#if defined(_DEBUG)
#  include <fmt/format.h>
#  include <fmt/ranges.h>
template <>
struct fmt::formatter<algo::edge<algo::detail::monostate>> {
  constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
  template <typename FormatContext>
  auto format(const algo::edge<algo::detail::monostate>& e,
              FormatContext& ctx) const {
    return fmt::format_to(ctx.out(), "(u: {}, v: {})", e.u, e.v);
  }
};
template <typename Edge>
struct fmt::formatter<algo::adjacency_list<Edge>> {
  constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
  template <typename FormatContext>
  auto format(const algo::adjacency_list<Edge>& adj, FormatContext& ctx) {
    auto buf = fmt::memory_buffer();
    fmt::format_to(back_inserter(buf), "(\n");
    for (int u = 0; u < adj.num_nodes(); ++u) {
      fmt::format_to(back_inserter(buf), "{}: [", u);
      for (auto it = adj[u].begin(); it != adj[u].end(); ++it) {
        fmt::format_to(back_inserter(buf), "(v: {}, w: {})", it->v, it->w);
        if (it != std::prev(adj[u].end())) {
          fmt::format_to(back_inserter(buf), ", ");
        }
      }
      fmt::format_to(back_inserter(buf), "]\n");
    }
    fmt::format_to(back_inserter(buf), ")");
    return fmt::format_to(ctx.out(), "{}", fmt::to_string(buf));
  }
  template <typename FormatContext>
  auto format(const algo::adjacency_list<Edge>& adj, FormatContext& ctx)
    requires std::same_as<algo::detail::monostate, Edge>
  {
    auto buf = fmt::memory_buffer();
    fmt::format_to(back_inserter(buf), "(\n");
    for (int u = 0; u < adj.num_nodes(); ++u) {
      fmt::format_to(back_inserter(buf), "{}: [", u);
      for (auto it = adj[u].begin(); it != adj[u].end(); ++it) {
        fmt::format_to(back_inserter(buf), "{}", *it);
        if (it != std::prev(adj[u].end())) {
          fmt::format_to(back_inserter(buf), ", ");
        }
      }
      fmt::format_to(back_inserter(buf), "]\n");
    }
    fmt::format_to(back_inserter(buf), ")");
    return fmt::format_to(ctx.out(), "{}", fmt::to_string(buf));
  }
};
#endif
