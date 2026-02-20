#pragma once

#include "../class/lazy-splay/point.hpp"

/*
@class/lazy-splay/point.hpp
*/

// makecode
namespace algo::lazy_splay {
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
auto make_point_query_handler(int n, Map map, Comp comp, CompId comp_id)
    -> point_query_handler<Tp, Mact, Map, Comp, CompId> {
  return {n, std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
auto make_point_query_handler(int n, const Tp& def, Map map, Comp comp,
                              CompId comp_id)
    -> point_query_handler<Tp, Mact, Map, Comp, CompId> {
  return {n, def, std::move(map), std::move(comp), std::move(comp_id)};
}
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId, detail::input_iterator<Tp> IIter>
auto make_point_query_handler(IIter first, IIter last, Map map, Comp comp,
                              CompId comp_id)
    -> point_query_handler<Tp, Mact, Map, Comp, CompId> {
  return {first, last, std::move(map), std::move(comp), std::move(comp_id)};
}
}  // namespace algo::lazy_splay