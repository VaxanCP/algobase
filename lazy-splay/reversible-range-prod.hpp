#pragma once

#include "../class/lazy-splay/reversible-range.hpp"

/*
@class/lazy-splay/reversible-range.hpp
*/
// makecode
namespace algo::lazy_splay {
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId>
auto make_reversible_range_query_handler(int n, Map map, Comp comp,
                                         CompId comp_id, BinOp binop,
                                         BinOpId binop_id)
    -> reversible_range_query_handler<Monoid, Mact, Map, Comp, CompId, BinOp,
                                      BinOpId> {
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
auto make_reversible_range_query_handler(int n, const Monoid& def, Map map,
                                         Comp comp, CompId comp_id, BinOp binop,
                                         BinOpId binop_id)
    -> reversible_range_query_handler<Monoid, Mact, Map, Comp, CompId, BinOp,
                                      BinOpId> {
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
auto make_reversible_range_query_handler(IIter first, IIter last, Map map,
                                         Comp comp, CompId comp_id, BinOp binop,
                                         BinOpId binop_id)
    -> reversible_range_query_handler<Monoid, Mact, Map, Comp, CompId, BinOp,
                                      BinOpId> {
  return {first,
          last,
          std::move(map),
          std::move(comp),
          std::move(comp_id),
          std::move(binop),
          std::move(binop_id)};
}
}  // namespace algo::lazy_splay
