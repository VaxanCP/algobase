#pragma once
#include <cassert>
#include <numeric>
// module header
#include "internal/base/def.hpp"
#include "internal/base/math-base.hpp"
/*
@internal/base/def.hpp
@internal/base/math-base.hpp
*/

// makecode
namespace algo {
/**
 * \brief Return floored square root of n
 *
 * \tparam Tp
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp sqrt(Tp n) {
#if !defined(NDEBUG)
  assert(n >= 0);
#endif
  return detail::sqrt(n);
}
/**
 * \brief Return absolute value of n.
 *
 * \tparam Tp
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp abs(Tp n) {
  return detail::abs(n);
}
/**
 * \brief Return floor(x / y). y must be non zero integer
 *
 * \tparam T1
 * \tparam T2
 * \param x
 * \param y y != 0 required
 * \return constexpr std::common_type_t<T1, T2>
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2>
constexpr auto floor_div(T1 x, T2 y) -> std::common_type_t<T1, T2> {
#if !defined(NDEBUG)
  assert(y != 0);
#endif
  using Tp = std::common_type_t<T1, T2>;
  return detail::floor_div<Tp>(x, y);
}
/**
 * \brief Return ceil(x / y). y must be non zero integer
 *
 * \tparam T1
 * \tparam T2
 * \param x
 * \param y
 * \return constexpr std::common_type_t<T1, T2>
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2>
constexpr auto ceil_div(T1 x, T2 y) -> std::common_type_t<T1, T2> {
#if !defined(NDEBUG)
  assert(y != 0);
#endif
  using Tp = std::common_type_t<T1, T2>;
  return detail::ceil_div<Tp>(x, y);
}
/**
 * \brief Return (-1)^n
 *
 * \tparam Tp
 * \param n
 * \return constexpr std::make_signed_t<std::common_type_t<Tp, char>>
 */
template <std::integral Tp>
constexpr int alt(Tp n) {
  return detail::alt(n);
}
// template <detail::integer T1, detail::sign_compatible_with<T1> T2,
//           detail::sign_compatible_with<T2> T3,
//           detail::sign_compatible_with<T3> T4>
// constexpr auto floor_sum(T1 a, T2 b, T3 c, T4 n)
//     -> std::common_type_t<T1, T2, T3, T4> {
// #if !defined(NDEBUG)
//   assert(n >= 0);
//   assert(c > 0);
//   assert(a >= 0);
// #endif
//   detail::assume(n >= 0);
//   detail::assume(c > 0);
//   detail::assume(a >= 0);
//   using Tp = std::common_type_t<T1, T2, T3, T4>;
//   Tp res = 0, a0 = a, b0 = b, c0 = c, n0 = n;
//   bool neg = false;
//   while (a0 != 0) {
//     if (a0 < c0 && b0 < c0) {
//       const Tp m = (a0 * n0 + b0) / c0;
//       const Tp tmp = a0;
//       res += alt(neg) * m * n0;
//       b0 = c0 - b0 - 1, a0 = c0, c0 = tmp, n0 = m - 1;
//       neg ^= true;
//     } else {
//       detail::assume(n0 >= 0);
//       const Tp tmp = (n0 * (n0 + 1) / 2) * (a0 / c0) + (n0 + 1) * (b0 / c0);
//       res += alt(neg) * tmp;
//       a0 %= c0, b0 %= c0;
//     }
//   }
//   res += alt(neg) * (b0 / c0) * (n0 + 1);
//   return res;
// }
}  // namespace algo