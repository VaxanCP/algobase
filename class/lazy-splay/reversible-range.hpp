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
template <typename Monoid, typename Mact>
struct reversible_range_query_node {
  using key_type = Monoid;
  reversible_range_query_node *children[2]{0};
  reversible_range_query_node *parent{0};
  Monoid self{};
  Monoid prod{};
  Mact lazy{};
  int size{1};
  bool prop{0};
  bool flip{0};
};
}  // namespace impl
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId>
class reversible_range_query_handler
    : public detail::basic_splay_tree<
          reversible_range_query_handler<Monoid, Mact, Map, Comp, CompId, BinOp,
                                         BinOpId>,
          impl::reversible_range_query_node<Monoid, Mact>> {
public:
  using self_type = reversible_range_query_handler;
  using node_type = impl::reversible_range_query_node<Monoid, Mact>;
  using base_type = detail::basic_splay_tree<self_type, node_type>;
  using typename base_type::const_node_pointer;
  using typename base_type::node_pointer;
  friend class detail::basic_splay_tree<self_type, node_type>;
  reversible_range_query_handler() noexcept = default;
  reversible_range_query_handler(int n, Map map, Comp comp, CompId comp_id,
                                 BinOp binop, BinOpId binop_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)),
        binop_(std::move(binop)),
        binop_id_(std::move(binop_id)) {
    this->root_ = build_tree(n, binop_id_());
  }
  reversible_range_query_handler(int n, const Monoid &def, Map map, Comp comp,
                                 CompId comp_id, BinOp binop,
                                 BinOpId binop_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)),
        binop_(std::move(binop)),
        binop_id_(std::move(binop_id)) {
    this->root_ = build_tree(n, def);
  }
  template <detail::input_iterator<Monoid> IIter>
  reversible_range_query_handler(IIter first, IIter last, Map map, Comp comp,
                                 CompId comp_id, BinOp binop,
                                 BinOpId binop_id) noexcept
      : base_type{},
        map_(std::move(map)),
        comp_(std::move(comp)),
        comp_id_(std::move(comp_id)),
        binop_(std::move(binop)),
        binop_id_(std::move(binop_id)) {
    const int n = static_cast<int>(std::distance(first, last));
    this->root_ = build_tree(first, n);
  }
  void apply(int i, const Mact &arg) {
#if !defined(NDEBUG)
    assert(0 <= i && i < this->size());
#endif
    apply(i, i + 1, arg);
  }
  void apply(int l, int r, const Mact &arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= this->size());
#endif
    if (l < r) {
      const auto node = this->access_range(l, r);
      apply_subtree(node, arg);
      this->reroot(node);
    }
  }
  void reverse(int l, int r) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= this->size());
#endif
    if (l < r) {
      const auto node = this->access_range(l, r);
      node->flip ^= 1;
      this->flip_range(node);
      this->reroot(node);
    }
  }
  Monoid query(int l, int r) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= this->size());
#endif
    // Ensure that call to access_range(l,r) is not a bug.
    if (l < r) [[likely]] {
      const auto node = this->access_range(l, r);
      const Monoid res = node->prod;
      this->reroot(node);
      return res;
    }
    return binop_id_();
  }
  Monoid query() const {
    if (!this->empty()) [[likely]] { return this->root_->prod; }
    return binop_id_();
  }
private:
  node_pointer build_tree(int n, const Monoid &def) const noexcept {
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
  void push_down(node_pointer u) const {
    flip_range(u);
    push_all(u);
  }
  /**
   * \brief
   *
   * \param u != nullptr
   */
  void push_up(node_pointer u) const {
    const auto [l, r] = u->children;
    Monoid prod_l = l != this->null() ? l->prod : binop_id_();
    const int size_l = this->size(l);
    Monoid prod_r = r != this->null() ? r->prod : binop_id_();
    const int size_r = this->size(r);
    Monoid sub_u = binop_(std::move(prod_l), std::move(prod_r));
    u->prod = binop_(u->self, std::move(sub_u));
    u->size = size_l + size_r + 1;
  }
  void flip_range(node_pointer u) const {
    if (u->flip) {
      std::swap(u->children[0], u->children[1]);
      const auto [l, r] = u->children;
      if (l != this->null()) { l->flip ^= 1; }
      if (r != this->null()) { r->flip ^= 1; }
      u->flip = false;
    }
  }
  void push_all(node_pointer u) const {
    if (u->prop) {
      const auto [l, r] = u->children;
      if (l != this->null()) { apply_subtree(l, u->lazy); }
      if (r != this->null()) { apply_subtree(r, u->lazy); }
      u->prop = false;
      u->lazy = comp_id_();
    }
  }
  // u is not null
  void apply_subtree(node_pointer u, const Mact &arg) const {
    u->self = map_(u->self, arg);
    u->prod = map_(u->prod, arg);
    u->lazy = comp_(u->lazy, arg);
    u->prop = true;
  }
  node_type make(Monoid key) const {
    return {.self = key, .prod = key, .lazy = comp_id_()};
  }
  [[no_unique_address]] Map map_;
  [[no_unique_address]] Comp comp_;
  [[no_unique_address]] CompId comp_id_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] BinOpId binop_id_;
};
}  // namespace algo::lazy_splay