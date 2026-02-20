#pragma once

#include "./class/lazy-segtree.hpp"

/*
@class/lazy-segtree.hpp
*/

// makecode
namespace algo {
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId>
auto make_lazy_segtree(int n, Map map, Comp comp, CompId comp_id, BinOp binop,
                       BinOpId binop_id)
    -> lazy_segtree<Monoid, Mact, Map, Comp, CompId, BinOp, BinOpId> {
  return {n,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(binop),
          std::move(binop_id)};
}
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId>
auto make_lazy_segtree(int n, const Monoid& def, Map map, Comp comp,
                       CompId comp_id, BinOp binop, BinOpId binop_id)
    -> lazy_segtree<Monoid, Mact, Map, Comp, CompId, BinOp, BinOpId> {
  return {n,
          def,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(binop),
          std::move(binop_id)};
}
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId,
          detail::input_iterator<Monoid> IIter>
auto make_lazy_segtree(IIter first, IIter last, Map map, Comp comp,
                       CompId comp_id, BinOp binop, BinOpId binop_id)
    -> lazy_segtree<Monoid, Mact, Map, Comp, CompId, BinOp, BinOpId> {
  return {first,
          last,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(binop),
          std::move(binop_id)};
}
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId, std::input_iterator IIter,
          detail::function<Monoid(std::iter_value_t<IIter>)> Gen>
auto make_lazy_segtree(IIter first, IIter last, Map map, Comp comp,
                       CompId comp_id, BinOp binop, BinOpId binop_id, Gen gen)
    -> lazy_segtree<Monoid, Mact, Map, Comp, CompId, BinOp, BinOpId> {
  return {first,
          last,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(binop),
          std::move(binop_id),
          std::move(gen)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
auto make_lazy_segtree(int n, Map map, Comp comp, CompId comp_id)
    -> lazy_segtree<Tp, Mact, Map, Comp, CompId> {
  return {n, std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
auto make_lazy_segtree(int n, const Tp& def, Map map, Comp comp, CompId comp_id)
    -> lazy_segtree<Tp, Mact, Map, Comp, CompId> {
  return {n, def, std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
auto make_lazy_segtree(std::vector<Tp> vec, Map map, Comp comp, CompId comp_id)
    -> lazy_segtree<Tp, Mact, Map, Comp, CompId> {
  return {std::move(vec), std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId, detail::input_iterator<Tp> IIter>
auto make_lazy_segtree(IIter first, IIter last, Map map, Comp comp,
                       CompId comp_id)
    -> lazy_segtree<Tp, Mact, Map, Comp, CompId> {
  return {first, last, std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId, std::input_iterator IIter,
          detail::function<Tp(std::iter_value_t<IIter>)> Gen>
auto make_lazy_segtree(IIter first, IIter last, Map map, Comp comp,
                       CompId comp_id, Gen gen)
    -> lazy_segtree<Tp, Mact, Map, Comp, CompId> {
  return {
      first,         last, std::move(map), std::move(comp), std::move(comp_id),
      std::move(gen)};
}
}  // namespace algo