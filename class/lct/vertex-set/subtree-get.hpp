#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "../../../internal/base/typing.hpp"
#include "../../../internal/lct.hpp"

/*
@internal/lct.hpp
@internal/base/typing.hpp
*/

// makecode
namespace algo::link_cut::vertex_set {
namespace impl {
template <typename Group>
struct subtree_aggregation_node {
  Group sub;
  Group self;
  Group vsub;
  int children[2]{0};
  int parent{0};
  bool flip{0};
};
}  // namespace impl
/**
 * \brief Link Cut Tree implementation for subtree aggregation queries.
 *
 * \tparam Group Abelian Group
 * \tparam BinOp Commutative Binary Operation
 * \tparam Id Identity Element
 * \tparam Inv Inverse Element
 */
template <std::semiregular Group, detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> Id, detail::function<Group(Group)> Inv>
class subtree_aggregation_handler
    : public detail::basic_lct<
          subtree_aggregation_handler<Group, BinOp, Id, Inv>,
          impl::subtree_aggregation_node<Group>> {
  using self_type = subtree_aggregation_handler;
  using node_type = impl::subtree_aggregation_node<Group>;
  using base_type = detail::basic_lct<self_type, node_type>;
  friend class detail::basic_lct<self_type, node_type>;
public:
  subtree_aggregation_handler() noexcept = default;
  subtree_aggregation_handler(int n, BinOp binop, Id id, Inv inv) noexcept
      : base_type{n, {id(), id(), id()}},
        binop_{std::move(binop)},
        id_{std::move(id)},
        inv_{std::move(inv)} {}
  subtree_aggregation_handler(int n, const Group& def, BinOp binop, Id id,
                              Inv inv) noexcept
      : base_type{n, {def, def, id()}},
        binop_{std::move(binop)},
        id_{std::move(id)},
        inv_{std::move(inv)} {
    tree_[Null].sub = id_();
    tree_[Null].self = id_();
  }
  template <detail::input_iterator<Group> IIter>
  subtree_aggregation_handler(IIter first, IIter last, BinOp binop, Id id,
                              Inv inv) noexcept
      : base_type{static_cast<int>(std::distance(first, last))},
        binop_{std::move(binop)},
        id_{std::move(id)},
        inv_{std::move(inv)} {
    std::transform(first, last, tree_.begin() + 1, [this](const Group& x) {
      return node_type{x, x, id_()};
    });
    tree_[Null].sub = id_();
    tree_[Null].self = id_();
    tree_[Null].vsub = id_();
  }
  /**
   * \brief Compute subtree-aggregate rooted at u
   *
   * \param u
   * \return
   */
  Group query_subtree(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    return get_subtree(u + 1);
  }
  /**
   * \brief Set new_val
   *
   * \param u
   * \param new_val
   */
  void set(int u, Group new_val) {
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
  Group get(int u) const {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    return do_get(u + 1);
  }
private:
  Group get_subtree(int u) {
    this->access(u);
    return binop_(tree_[u].self, tree_[u].vsub);
  }
  void do_set(int u, Group new_val) {
    this->access(u);
    tree_[u].self = std::move(new_val);
    this->update(u);
  }
  Group do_get(int u) const { return tree_[u].self; }
  void push_down(int u) {
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
    Group sub = binop_(tree_[l].sub, tree_[r].sub);
    tree_[u].sub = binop_(binop_(tree_[u].vsub, tree_[u].self), std::move(sub));
  }
  void fix_remove(int v, int c) {
    tree_[v].vsub = binop_(tree_[v].vsub, tree_[c].sub);
  }
  void fix_attach(int v, int u) {
    tree_[v].vsub = binop_(tree_[v].vsub, inv_(tree_[u].sub));
  }
  using base_type::Null;
  using base_type::tree_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
  [[no_unique_address]] Inv inv_;
};
}  // namespace algo::link_cut::vertex_set