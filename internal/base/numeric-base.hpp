#pragma once
#include <array>
#include <numeric>
#include <utility>

#include "./def.hpp"
#include "./math-base.hpp"
#include "./typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/def.hpp
@internal/base/math-base.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Safe integer multiplication
 *
 * \tparam Tp integer required
 * \param x
 * \param y
 * \return constexpr imul_result_t<Tp>
 */
template <typename Tp>
constexpr auto safe_mul(Tp x, Tp y) -> imul_result_t<Tp> {
  using promoted_type = imul_result_t<Tp>;
  return static_cast<promoted_type>(x) * y;
}
/**
 * \brief Compute x mod m.
 *
 * \tparam Tp
 * \param x
 * \param m 0 < m required.
 * \return constexpr Tp Integer n with 0 <= n < m such that n == x mod m.
 */
template <typename Tp>
constexpr Tp norm_mod(Tp x, Tp m) {
  if ((x %= m) < 0) { x += m; }
  return x;
}
/**
 * \brief modular multiplication
 *
 * \tparam Tp Integer required
 * \param x 0 <= x < m required
 * \param y 0 <= y < m required
 * \param m Modulus. 0 < m required
 * \return constexpr Tp Integer z with z = x*y mod m.
 */
template <typename Tp>
constexpr Tp mul_mod(Tp x, Tp y, Tp m) {
  assume(0 < m);
  assume(0 <= x && x < m);
  assume(0 <= y && y < m);
  const auto prod = safe_mul(x, y) % m;
  return static_cast<Tp>(prod);
}
/**
 * \brief Extended euclid algorithm. Return any solution (x,y) of the form
 * ax+by = gcd(a,b)
 *
 * \tparam Tp Signed integer required
 * \param a
 * \param b
 * \return
 */
template <signed_integer Tp>
constexpr auto euclid_algo(Tp a, Tp b) -> std::array<Tp, 3> {
  Tp x0 = 1, y0 = 0;
  Tp x1 = 0, y1 = 1;
  Tp a0 = a, b0 = b;
  while (b0) {
    Tp q = a0 / b0, tp = x0;
    x0 = x1;
    x1 = tp - q * x1;
    tp = y0;
    y0 = y1;
    y1 = tp - q * y1;
    tp = a0;
    a0 = b0;
    b0 = tp - q * b0;
  }
  if (a0 < 0) { x0 = -x0, y0 = -y0, a0 = -a0; }
  return {x0, y0, a0};
}
/**
 * \brief Perform binary exponentiation.
 *
 * \tparam Tp
 * \tparam Up Integer required.
 * \param b Base.
 * \param n Exponent. 0 <= n required.
 * \param id Multiplicative identity.
 * \param binop Binary Operation
 * \return constexpr Tp
 */
template <typename Tp, typename Up, typename BinOp = std::multiplies<Tp>>
constexpr Tp binary_exp(Tp b, Up n, std::type_identity_t<Tp> id,
                        BinOp binop = {}) {
  Tp ans = id;
  while (n > 0) {
    if (n % 2 == 1) { ans = binop(ans, b); }
    b = binop(b, b);
    n /= 2;
  }
  return ans;
}
/**
 * \brief Compute modular exponentiation.
 *
 * \tparam Tp Integer required.
 * \tparam Up Integer required.
 * \param b Base. 0 <= b < m required.
 * \param n Exponent. 0 <= n required.
 * \param m Modulus. m > 0 required.
 * \return constexpr Tp b^n mod m.
 */
template <typename Tp, typename Up>
constexpr Tp pow_mod(Tp b, Up n, Tp m) {
  assume(m > 0);
  if (m == 1) [[unlikely]] { return 0; }
  return binary_exp(b, n, 1, [m](Tp x, Tp y) { return mul_mod(x, y, m); });
}
}  // namespace algo::detail