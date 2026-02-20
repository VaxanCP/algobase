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
namespace algo::link_cut::subtree_set {
namespace impl {
template <typename Group, typename Gact>
struct path_subtree_aggregation_node {
  Group self;  // value at a node
  Group sub;   // subtree aggregation
  Group path;  // path aggregation
  Group vsub;  // virtual subtree aggregation
  Gact lazy;   // lazy value
  int children[2]{0};
  int parent{0};
  bool flip{0};
};
}  // namespace impl
/**
 * \brief Link Cut Tree implementation for subtree update and path and subtree
 * aggregation queries.
 *
 * \tparam Group
 * \tparam Gact
 * \tparam Map
 * \tparam Comp
 * \tparam CompId
 * \tparam CompInv
 * \tparam BinOp
 * \tparam BinOpId
 * \tparam BinOpInv
 */
template <std::semiregular Group, std::semiregular Gact,
          detail::function<Group(Group, Gact)> Map,
          detail::function<Gact(Gact, Gact)> Comp,
          detail::function<Gact(void)> CompId,
          detail::function<Gact(Gact)> CompInv,
          detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> BinOpId,
          detail::function<Group(Group)> BinOpInv>
class path_subtree_aggregation_handler
    : public detail::basic_lct<
          path_subtree_aggregation_handler<Group, Gact, Map, Comp, CompId,
                                           CompInv, BinOp, BinOpId, BinOpInv>,
          impl::path_subtree_aggregation_node<Group, Gact>> {
  using self_type = path_subtree_aggregation_handler;
  using node_type = impl::path_subtree_aggregation_node<Group, Gact>;
  using base_type = detail::basic_lct<self_type, node_type>;
  friend class detail::basic_lct<self_type, node_type>;
public:
  path_subtree_aggregation_handler() noexcept = default;
  path_subtree_aggregation_handler(int n, Map map, Comp comp, CompId comp_id,
                                   CompInv comp_inv, BinOp binop,
                                   BinOpId binop_id,
                                   BinOpInv binop_inv) noexcept
      : base_type{n,
                  {binop_id(), binop_id(), binop_id(), binop_id(), comp_id()}},
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)},
        comp_inv_{std::move(comp_inv)},
        binop_{std::move(binop)},
        binop_id_{std::move(binop_id)},
        binop_inv_{std::move(binop_inv)} {}
  path_subtree_aggregation_handler(int n, const Group& def, Map map, Comp comp,
                                   CompId comp_id, CompInv comp_inv,
                                   BinOp binop, BinOpId binop_id,
                                   BinOpInv binop_inv) noexcept
      : base_type{n, {def, def, def, binop_id(), comp_id()}},
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)},
        comp_inv_{std::move(comp_inv)},
        binop_{std::move(binop)},
        binop_id_{std::move(binop_id)},
        binop_inv_{std::move(binop_inv)} {
    tree_[Null].self = binop_id_();
    tree_[Null].sub = binop_id_();
    tree_[Null].path = binop_id_();
  }
  Group query_subtree(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    return get_subtree(u + 1);
  }
  Group query_path(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
    assert(0 <= v && v < this->num_nodes());
    assert(this->connected(u, v));
#endif
    return get_path(v + 1, u + 1);
  }
  void apply_subtree(int u, Gact z) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    act_subtree(u + 1, std::move(z));
  }
  void set(int u, Group new_val) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    do_set(u + 1, std::move(new_val));
  }
  Group get(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < this->num_nodes());
#endif
    return do_get(u + 1);
  }
private:
  Group get_subtree(int u) {
    this->access(u);
    Group sub_u = binop_(tree_[u].vsub, tree_[u].self);
    return map_(std::move(sub_u), tree_[u].lazy);
  }
  Group get_path(int v, int u) {
    this->do_reroot(u);
    this->access(v);
    return map_(tree_[v].path, tree_[v].lazy);
  }
  void act_subtree(int u, Gact z) {
    // Can path update be possible with non-commutative action?
    this->access(u);
    tree_[u].lazy = comp_(tree_[u].lazy, z);
    const int l = tree_[u].children[0];
    if (l != Null) {
      tree_[l].lazy = comp_(tree_[l].lazy, comp_inv_(std::move(z)));
      this->update(u);
    }
  }
  void do_set(int u, Group new_val) {
    // Normal point update
    this->access(u);
    tree_[u].self = std::move(new_val);
    this->update(u);
  }
  Group do_get(int u) {
    this->access(u);
    return map_(tree_[u].self, tree_[u].lazy);
  }
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
    Group sub_l = map_(tree_[l].sub, tree_[l].lazy);
    Group sub_r = map_(tree_[r].sub, tree_[r].lazy);
    Group path_l = map_(tree_[l].path, tree_[l].lazy);
    Group path_r = map_(tree_[r].path, tree_[r].lazy);
    tree_[u].sub = binop_(binop_(tree_[u].self, tree_[u].vsub),
                          binop_(std::move(sub_l), std::move(sub_r)));
    tree_[u].path =
        binop_(tree_[u].self, binop_(std::move(path_l), std::move(path_r)));
  }
  void fix_link(int v, int u) {
    tree_[u].lazy = comp_(tree_[u].lazy, comp_inv_(tree_[v].lazy));
  }
  void fix_cut(int v, int u) {
    tree_[u].lazy = comp_(tree_[u].lazy, tree_[v].lazy);
  }
  void fix_rotate(int v, int u, int b) {
    // d`[u] = d[u] * d[v]
    // d`[v] = d[u]^(-1)
    // d`[b] = d[b] * d[u]
    // Refer to:
    // http://www.planarity.org/Klein_splay_trees_and_link-cut_trees.pdf
    tree_[u].lazy = comp_(tree_[u].lazy, tree_[v].lazy);
    tree_[v].lazy = comp_(tree_[v].lazy, comp_inv_(tree_[u].lazy));
    if (b != Null) {
      tree_[b].lazy = comp_(tree_[b].lazy, comp_inv_(tree_[v].lazy));
    }
  }
  void fix_remove(int v, int c) {
    // Refer to:
    // http://www.planarity.org/Klein_splay_trees_and_link-cut_trees.pdf
    Group sub_c = map_(tree_[c].sub, tree_[c].lazy);
    tree_[v].vsub = binop_(tree_[v].vsub, std::move(sub_c));
  }
  void fix_attach(int v, int u) {
    // Refer to:
    // http://www.planarity.org/Klein_splay_trees_and_link-cut_trees.pdf
    Group sub_u = map_(tree_[u].sub, tree_[u].lazy);
    tree_[v].vsub = binop_(tree_[v].vsub, binop_inv_(std::move(sub_u)));
  }
  using base_type::Null;
  using base_type::tree_;
  [[no_unique_address]] Map map_;
  [[no_unique_address]] Comp comp_;
  [[no_unique_address]] CompId comp_id_;
  [[no_unique_address]] CompInv comp_inv_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] BinOpId binop_id_;
  [[no_unique_address]] BinOpInv binop_inv_;
};
}  // namespace algo::link_cut::subtree_set