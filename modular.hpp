#pragma once
#include <cassert>
#include <numeric>
// modules
#include "internal/base/def.hpp"
#include "internal/base/typing.hpp"
#include "internal/modular.hpp"
#include "internal/numeric.hpp"
/*
@internal/base/typing.hpp
@internal/base/def.hpp
@internal/numeric.hpp
@internal/modular.hpp
*/

// makecode
namespace algo {
/**
 * \brief Compute modular exponentiation.
 *
 * \tparam T1
 * \tparam T2
 * \tparam Up Integer required.
 * \param b Base. 0 <= b < m required.
 * \param n Exponent. 0 <= n required.
 * \param m Modulus. 0 < m required.
 * \return std::common_type_t<T1, T2> b^n.
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2,
          detail::integer Up>
constexpr auto pow_mod(T1 b, Up n, T2 m) -> std::common_type_t<T1, T2> {
#if !defined(NDEBUG)
  assert(m > 0);
  assert(0 <= b && b < m);
  assert(n >= 0);
#endif
  using Tp = std::common_type_t<T1, T2>;
  return detail::pow_mod<Tp, Up>(b, n, m);
}
/**
 * \brief Return modular multiplicative inverse
 *
 * \tparam Tp
 * \param a
 * \param m m > 0 required.
 * \return constexpr Tp
 * \note gcd(a,m) = 1 should hold
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2>
constexpr auto inv_mod(T1 a, T2 m) -> std::common_type_t<T1, T2> {
#if !defined(NDEBUG)
  assert(std::gcd(a, m) == 1);
  assert(m > 0);
#endif
  using Tp = std::common_type_t<T1, T2>;
  return detail::inv_mod<Tp>(a, m);
}
/**
 * \brief Find a primitive root modulo n or nullopt if it doesn't exist.
 Primitive root modulo n exists if and only if:
 1. n is 1,2,4 or
 2. n is power of an odd prime number (i.e. n=p^k for some prime number p)
 or
 3. n is of the form 2p^k for some odd prime number p
 time complexity: O(Ans*log(phi(n))*log(n)) where phi(x) is euler totient
 function
 *
 * \tparam Tp
 * \param n
 * \param phi The number of positive integer coprime to n (i.e. The value
 of euler totient function at n)
 * \return constexpr std::optional<Tp>
 */
template <detail::integer Tp>
constexpr auto primitive_root(Tp n, Tp phi) -> std::optional<Tp> {
#if !defined(NDEBUG)
  assert(n > 0);
  assert(phi == detail::totient(n));
#endif
  return detail::primitive_root(n, phi);
}
/**
 * \brief Return primitive root modulo n where n is prime
 *
 * \tparam Tp
 * \param n prime number
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp primitive_root_prime(Tp n) {
#if !defined(NDEBUG)
  assert(detail::prime_test(n));
#endif
  return detail::primitive_root_prime(n);
}
/**
 * \brief find one solution for a^x == b (mod n) using Baby-step giant-step
algorithm time complexity: O(sqrt(n))
 *
 * \tparam T1
 * \tparam T2
 * \tparam T3
 * \param a 0 <= a < n required
 * \param b 0 <= b < n required
 * \param n n > 0 required.
 * \return std::optional<std::common_type_t<T1, T2, T3>>
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2,
          detail::sign_compatible_with<T2> T3>
auto discrete_log(T1 a, T2 b, T3 n)
    -> std::optional<std::common_type_t<T1, T2, T3>> {
#if !defined(NDEBUG)
  assert(n > 0);
  assert(0 <= a && a < n);
  assert(0 <= b && b < n);
#endif
  using Tp = std::common_type_t<T1, T2, T3>;
  return detail::discrete_log<Tp>(a, b, n);
}
/**
 * \brief find one solution for x: x^a == b mod n (where n is prime)
if a == 0 and b == 0, result is undefined
for general solution, solutions can be represented as follows:
ans=g^(y0+(i*(phi(n)/gcd(k, phi(n))))) mod n where i is some integer and g
is primitive root modulo n. time complexity: O(sqrt(n)+g*log(n)^2)
 *
 * \tparam T1
 * \tparam T2
 * \tparam T3
 * \param a
 * \param b 0 <= b < n required.
 * \param n prime number
 * \return std::optional<std::common_type_t<T1, T2, T3>>
 */
template <detail::integer T1, detail::sign_compatible_with<T1> T2,
          detail::sign_compatible_with<T2> T3>
auto discrete_root(T1 a, T2 b, T3 n)
    -> std::optional<std::common_type_t<T1, T2, T3>> {
#if !defined(NDEBUG)
  assert(detail::prime_test(n));
  assert(a >= 0);
  assert(0 <= b && b < n);
#endif
  using Tp = std::common_type_t<T1, T2, T3>;
  return detail::discrete_root<Tp>(a, b, n);
}
}  // namespace algo