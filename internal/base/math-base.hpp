#pragma once
#include <bit>
#include <numeric>

#include "./bit-base.hpp"
#include "./typing.hpp"
/*
@internal/base/typing.hpp
@internal/base/bit-base.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Do what you think it should do
 *
 * \tparam Tp Integer required
 * \param n n >= 0 required
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp sqrt(Tp n) {
  if (n == 0) { return 0; }
  using Up = std::make_unsigned_t<Tp>;
  // To avoild overflow.
  const Up x = static_cast<Up>(n);
  Up ans = 0;
  for (int i = floor_log2(n) / 2; i >= 0; --i) {
    const Up tmp = ans | (Up(1) << i);
    if (tmp * tmp <= x) { ans = tmp; }
  }
  return static_cast<Tp>(ans);
}
/**
 * \brief Return absolute value of the given value.
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp abs(Tp n) {
  return n < 0 ? -n : n;
}
/**
 * \brief Return floor(x / y)
 *
 * \tparam Tp Integer required
 * \param x
 * \param y y != 0 required
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp floor_div(Tp x, Tp y) {
  if constexpr (std::is_signed_v<Tp>) {
    return x / y - ((x ^ y) < 0 && x % y != 0);
  } else {
    return x / y;
  }
}
/**
 * \brief Return ceil(x / y)
 *
 * \tparam Tp Integer required
 * \param x
 * \param y y != 0 required
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp ceil_div(Tp x, Tp y) {
  if constexpr (std::is_signed_v<Tp>) {
    return y > 0 ? floor_div(x + y - 1, y) : -floor_div(x, -y);
  } else {
    return floor_div(x + y - 1, y);
  }
}
/**
 * \brief Return (-1)^n
 *
 * \tparam Tp Integral required
 * \param n
 * \return constexpr std::make_signed_t<std::common_type_t<Tp, int>>
 */
template <typename Tp>
constexpr int alt(Tp n) {
  return -(n & 1) | 1;
}
}  // namespace algo::detail