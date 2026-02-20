#pragma once
#include <algorithm>
#include <cassert>
#include <vector>

#include "../../../internal/base/typing.hpp"
#include "../../../internal/lct.hpp"

/*
@internal/base/typing.hpp
@internal/lct.hpp
*/

// makecode
namespace algo::link_cut::vertex_set {
namespace impl {
template <typename Monoid>
struct path_aggregation_node {
  Monoid path[2];
  Monoid self;
  int children[2]{0};
  int parent{0};
  bool flip{0};
};
}  // namespace impl
/**
 * \brief Link Cut Tree implementation for point update and path query.
 *
 * \tparam Monoid Element type
 * \tparam BinOp Binary Operation. This doesn't need to be commutative.
 * \tparam Id Identity Element.
 */
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
class path_aggregation_handler
    : public detail::basic_lct<path_aggregation_handler<Monoid, BinOp, Id>,
                               impl::path_aggregation_node<Monoid>> {
  using self_type = path_aggregation_handler;
  using node_type = impl::path_aggregation_node<Monoid>;
  using base_type = detail::basic_lct<self_type, node_type>;
  friend class detail::basic_lct<self_type, node_type>;
public:
  path_aggregation_handler() noexcept = default;
  path_aggregation_handler(int n, BinOp binop, Id id) noexcept
      : base_type{n, {id(), id(), id()}},
        binop_{std::move(binop)},
        id_{std::move(id)} {}
  path_aggregation_handler(int n, const Monoid& def, BinOp binop,
                           Id id) noexcept
      : base_type{n, {def, def, def}},
        binop_{std::move(binop)},
        id_{std::move(id)} {
    tree_[Null].path[0] = id_();
    tree_[Null].path[1] = id_();
    tree_[Null].self = id_();
  }
  template <detail::input_iterator<Monoid> IIter>
  path_aggregation_handler(IIter first, IIter last, BinOp binop, Id id) noexcept
      : base_type{static_cast<int>(std::distance(first, last))},
        binop_{std::move(binop)},
        id_{std::move(id)} {
    std::transform(first, last, tree_.begin() + 1, [](const Monoid& x) {
      return node_type{x, x, x};
    });
    // Set the aggregation values of null node to identity.
    tree_[Null].path[0] = id_();
    tree_[Null].path[1] = id_();
    tree_[Null].self = id_();
  }
  /**
   * \brief Compute path-aggregate between u and v. u and v must be
   * connected.
   *
   * \param v
   * \param u
   * \return
   */
  Monoid query_path(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
    assert(0 <= v && v < this->num_nodes());
    assert(this->connected(u, v));
#endif
    return get_path(v + 1, u + 1);
  }
  /**
   * \brief Set new_val
   *
   * \param u
   * \param new_val
   */
  void set(int u, Monoid new_val) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    do_set(u + 1, std::move(new_val));
  }
  /**
   * \brief Get the value of u
   *
   * \param u
   * \return
   */
  Monoid get(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    return do_get(u + 1);
  }
private:
  Monoid get_path(int v, int u) {
    this->do_reroot(v);
    this->access(u);
    return tree_[u].path[0];
  }
  void do_set(int u, Monoid new_val) {
    this->access(u);
    tree_[u].self = std::move(new_val);
    this->update(u);
  }
  Monoid do_get(int u) const { return tree_[u].self; }
  void push_down(int u) {
    // If the node u is to be flipped, do flip the childrens and propagate the
    // tag to them.
    if (tree_[u].flip) {
      std::swap(tree_[u].children[0], tree_[u].children[1]);
      const int l = tree_[u].children[0];
      const int r = tree_[u].children[1];
      tree_[u].flip ^= true;
      tree_[l].flip ^= true;
      tree_[r].flip ^= true;
    }
  }
  void push_up(int u) {
    const int l = tree_[u].children[0];
    const int r = tree_[u].children[1];
    const bool fl = tree_[l].flip;
    const bool fr = tree_[r].flip;
    tree_[u].path[0] = binop_(tree_[l].path[fl], tree_[u].self);
    tree_[u].path[0] = binop_(tree_[u].path[0], tree_[r].path[fr]);
    tree_[u].path[1] = binop_(tree_[r].path[!fr], tree_[u].self);
    tree_[u].path[1] = binop_(tree_[u].path[1], tree_[l].path[!fl]);
  }
  using base_type::Null;
  using base_type::tree_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
};
}  // namespace algo::link_cut::vertex_set