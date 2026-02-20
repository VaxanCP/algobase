#pragma once

#include "class/sparse-table.hpp"

/*
@class/sparse-table.hpp
*/

// makecode
namespace algo {
template <std::semiregular SemiGroup,
          detail::function<SemiGroup(SemiGroup, SemiGroup)> BinOp>
auto make_sparse_table(int n, const SemiGroup& def, BinOp binop)
    -> sparse_table<SemiGroup, BinOp> {
  return {n, def, std::move(binop)};
}
template <std::semiregular SemiGroup,
          detail::function<SemiGroup(SemiGroup, SemiGroup)> BinOp,
          detail::input_iterator<SemiGroup> IIter>
auto make_sparse_table(IIter first, IIter last, BinOp binop)
    -> sparse_table<SemiGroup, BinOp> {
  return {first, last, std::move(binop)};
}
template <std::semiregular SemiGroup,
          detail::function<SemiGroup(SemiGroup, SemiGroup)> BinOp,
          std::input_iterator IIter,
          detail::function<SemiGroup(std::iter_value_t<IIter>)> Gen>
auto make_sparse_table(IIter first, IIter last, BinOp binop, Gen gen)
    -> sparse_table<SemiGroup, BinOp> {
  return {first, last, std::move(binop), std::move(gen)};
}
}  // namespace algo