#pragma once

#include "./class/persegtree.hpp"

/*
@class/persegtree.hpp
*/

// makecode
namespace algo {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_persistent_segtree(int n, BinOp binop, Id id)
    -> persistent_segtree<Monoid, BinOp, Id> {
  return {n, std::move(binop), std::move(id)};
}
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_persistent_segtree(int n, const Monoid& def, BinOp binop, Id id)
    -> persistent_segtree<Monoid, BinOp, Id> {
  return {n, def, std::move(binop), std::move(id)};
}
template <
    std::semiregular Monoid, detail::function<Monoid(Monoid, Monoid)> BinOp,
    detail::function<Monoid(void)> Id, detail::input_iterator<Monoid> IIter>
auto make_persistent_segtree(IIter first, IIter last, BinOp binop, Id id)
    -> persistent_segtree<Monoid, BinOp, Id> {
  return {first, last, std::move(binop), std::move(id)};
}
}  // namespace algo