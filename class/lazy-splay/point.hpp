#pragma once

#include <cassert>
#include <memory_resource>
#include <utility>

#include "../../internal/base/typing.hpp"
#include "../../internal/splay.hpp"

/*
@internal/splay.hpp
@internal/base/typing.hpp
*/

// makecode
namespace algo::lazy_splay {
namespace impl {
template <typename Tp, typename Mact>
struct point_query_node {
  using key_type = Tp;
  point_query_node* children[2]{0};
  point_query_node* parent{0};
  Tp self{};
  Mact lazy{};
  int size{1};
};
}  // namespace impl
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
class point_query_handler
    : public detail::basic_splay_tree<
          point_query_handler<Tp, Mact, Map, Comp, CompId>,
          impl::point_query_node<Tp, Mact>> {
public:
  using self_type = point_query_handler;
  using node_type = impl::point_query_node<Tp, Mact>;
  using base_type = detail::basic_splay_tree<self_type, node_type>;
  using typename base_type::const_node_pointer;
  using typename base_type::node_pointer;
  friend class detail::basic_splay_tree<self_type, node_type>;
  point_query_handler() noexcept = default;
  point_query_handler(int n, Map map, Comp comp, CompId comp_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)) {
    this->root_ = build_tree(n, Tp{});
  }
  point_query_handler(int n, const Tp& def, Map map, Comp comp,
                      CompId comp_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)) {
    this->root_ = build_tree(n, def);
  }
  template <detail::input_iterator<Tp> IIter>
  point_query_handler(IIter first, IIter last, Map map, Comp comp,
                      CompId comp_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)) {
    const int n = static_cast<int>(std::distance(first, last));
    this->root_ = build_tree(first, n);
  }
  void apply(int i, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= i && i < this->size());
#endif
    apply(i, i + 1, arg);
  }
  void apply(int l, int r, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= this->size());
#endif
    if (l < r) {
      const auto node = this->access_range(l, r);
      apply_subtree(node, arg);
      this->reroot(node);
    }
  }
protected:
  node_pointer build_tree(int n, const Tp& def) const noexcept {
    if (n == 0) {
      return this->null();
    } else {
      const auto node = this->allocate(def);
      this->attach(node, build_tree(n >> 1, def), 0);
      this->attach(node, build_tree((n - 1) >> 1, def), 1);
      return node;
    }
  }
  template <typename IIter>
  node_pointer build_tree(IIter base, int n) const noexcept {
    if (n == 0) {
      return this->null();
    } else {
      const auto mid = base + (n >> 1);
      const auto node = this->allocate(*mid);
      this->attach(node, build_tree(base, n >> 1), 0);
      this->attach(node, build_tree(mid + 1, (n - 1) >> 1), 1);
      return node;
    }
  }
  void push_up(node_pointer u) const {
    const auto [l, r] = u->children;
    const int size_l = this->size(l);
    const int size_r = this->size(r);
    u->size = size_l + size_r + 1;
  }
  // u != nullptr
  void push_down(node_pointer u) const {
    const auto [l, r] = u->children;
    if (l != this->null()) { apply_subtree(l, u->lazy); }
    if (r != this->null()) { apply_subtree(r, u->lazy); }
    u->lazy = comp_id_();
  }
  // u != nullptr
  void apply_subtree(node_pointer u, const Mact& arg) const {
    u->self = map_(u->self, arg);
    u->lazy = comp_(u->lazy, arg);
  }
  /**
   * \brief Make node based on key.
   *
   * \param key
   * \return
   */
  node_type make(Tp key) const {
    return {.self = std::move(key), .lazy = comp_id_()};
  }
  [[no_unique_address]] Map map_;
  [[no_unique_address]] Comp comp_;
  [[no_unique_address]] CompId comp_id_;
};
}  // namespace algo::lazy_splay