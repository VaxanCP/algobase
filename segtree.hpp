#pragma once

#include "./class/segtree.hpp"

/*
@class/segtree.hpp
*/

// makecode
namespace algo {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_segtree(int n, BinOp binop, Id id) -> segtree<Monoid, BinOp, Id> {
  return {n, std::move(binop), std::move(id)};
}
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_segtree(int n, const Monoid& def, BinOp binop, Id id)
    -> segtree<Monoid, BinOp, Id> {
  return {n, def, std::move(binop), std::move(id)};
}
template <
    std::semiregular Monoid, detail::function<Monoid(Monoid, Monoid)> BinOp,
    detail::function<Monoid(void)> Id, detail::input_iterator<Monoid> IIter>
auto make_segtree(IIter first, IIter last, BinOp binop, Id id)
    -> segtree<Monoid, BinOp, Id> {
  return {first, last, std::move(binop), std::move(id)};
}
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id, std::input_iterator IIter,
          detail::function<Monoid(std::iter_value_t<IIter>)> Gen>
auto make_segtree(IIter first, IIter last, BinOp binop, Id id, Gen gen)
    -> segtree<Monoid, BinOp, Id> {
  return {first, last, std::move(binop), std::move(id), std::move(gen)};
}
}  // namespace algo