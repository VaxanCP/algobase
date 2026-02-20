#pragma once
#include "class/fenwick-tree.hpp"

/*
@class/fenwick-tree.hpp
*/
// makecode
namespace algo {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_fenwick_tree(int n, BinOp binop, Id id)
    -> fenwick_tree<Monoid, BinOp, Id> {
  return {n, std::move(binop), std::move(id)};
}
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_fenwick_tree(int n, const Monoid& def, BinOp binop, Id id)
    -> fenwick_tree<Monoid, BinOp, Id> {
  return {n, def, std::move(binop), std::move(id)};
}
}  // namespace algo