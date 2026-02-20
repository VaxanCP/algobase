#pragma once
#include <cassert>
// makecode
namespace algo::detail {
/**
 * \brief
 *
 */
inline void unreachable [[gnu::noreturn]] () { __builtin_unreachable(); }
/**
 * \brief It can enable some optimizations based on the assumption that the
 * given expression can never evaluate to false. Note that if the expression is
 * false, the behavior is undefined. Does runtime check if NDEBUG is not
 * defined.
 *
 * \param expr
 */
constexpr void assume [[gnu::always_inline]] (bool expr) {
#if !defined(NDEBUG)
  assert(expr);
#endif
  if (!expr) { __builtin_unreachable(); }
}
}  // namespace algo::detail
