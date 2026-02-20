#pragma once

#include "./class/swag-queue.hpp"

/*
@class/swag-queue.hpp
*/
// makecode
namespace algo {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
auto make_swag_queue(BinOp binop, Id id) -> swag_queue<Monoid, BinOp, Id> {
  return {std::move(binop), std::move(id)};
}
}  // namespace algo