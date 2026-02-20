#pragma once

#include "../../class/lct/vertex-set/subtree-get.hpp"

/*
@class/lct/vertex-set/subtree-get.hpp
*/

// makecode
namespace algo::link_cut::vertex_set {
template <std::semiregular Group, detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> Id, detail::function<Group(Group)> Inv>
auto make_subtree_aggregation_handler(int n, BinOp binop, Id id, Inv inv)
    -> subtree_aggregation_handler<Group, BinOp, Id, Inv> {
  return {n, std::move(binop), std::move(id), std::move(inv)};
}
template <std::semiregular Group, detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> Id, detail::function<Group(Group)> Inv>
auto make_subtree_aggregation_handler(int n, const Group& def, BinOp binop,
                                      Id id, Inv inv)
    -> subtree_aggregation_handler<Group, BinOp, Id, Inv> {
  return {n, def, std::move(binop), std::move(id), std::move(inv)};
}
template <std::semiregular Group, detail::function<Group(Group, Group)> BinOp,
          detail::function<Group(void)> Id, detail::function<Group(Group)> Inv,
          detail::input_iterator<Group> IIter>
auto make_subtree_aggregation_handler(IIter first, IIter last, BinOp binop,
                                      Id id, Inv inv)
    -> subtree_aggregation_handler<Group, BinOp, Id, Inv> {
  return {first, last, std::move(binop), std::move(id), std::move(inv)};
}
}  // namespace algo::link_cut::vertex_set