#pragma once

#include "../../class/lct/subtree-set/path-subtree-get.hpp"

/*
@class/lct/subtree-set/path-subtree-get.hpp
*/

// makecode
namespace algo::link_cut::subtree_set {
template <std::semiregular Group, std::semiregular Gact,
          detail::function<Group(Group, Gact)> Map,
          detail::function<Gact(Gact, Gact)> Comp,
          detail::function<Gact(void)> CompId,
          detail::function<Gact(Gact)> CompInv,
          detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> BinOpId,
          detail::function<Group(Group)> BinOpInv>
auto make_path_subtree_aggregation_handler(int n, Map map, Comp comp,
                                           CompId comp_id, CompInv comp_inv,
                                           BinOp binop, BinOpId binop_id,
                                           BinOpInv binop_inv) {
  return {n,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(comp_inv),
          std::move(binop),
          std::move(binop_id),
          std::move(binop_inv)};
}
}  // namespace algo::link_cut::subtree_set