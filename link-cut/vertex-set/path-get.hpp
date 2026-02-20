#pragma once

#include "../../class/lct/vertex-set/path-get.hpp"

/*
@class/lct/vertex-set/path-get.hpp
*/

// makecode
namespace algo::link_cut::vertex_set {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_path_aggregation_handler(int n, BinOp binop, Id id)
    -> path_aggregation_handler<Monoid, BinOp, Id> {
  return {n, std::move(binop), std::move(id)};
}
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_path_aggregation_handler(int n, const Monoid& def, BinOp binop, Id id)
    -> path_aggregation_handler<Monoid, BinOp, Id> {
  return {n, def, std::move(binop), std::move(id)};
}
template <
    std::semiregular Monoid, detail::function<Monoid(Monoid, Monoid)> BinOp,
    detail::function<Monoid(void)> Id, detail::input_iterator<Monoid> IIter>
auto make_path_aggregation_handler(IIter first, IIter last, BinOp binop, Id id)
    -> path_aggregation_handler<Monoid, BinOp, Id> {
  return {first, last, std::move(binop), std::move(id)};
}
}  // namespace algo::link_cut::vertex_set