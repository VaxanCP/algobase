#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <optional>
// modules
#include "base/def.hpp"
#include "base/math-base.hpp"
#include "base/numeric-base.hpp"
#include "base/typing.hpp"
/*
@internal/base/math-base.hpp
@internal/base/numeric-base.hpp
@internal/base/typing.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Return modular multiplicative inverse
 * \tparam Tp Integer required
 * \param a gcd(a, m) = 1 required
 * \param m m > 0 required
 * \return constexpr Tp 0 <= res < m will hold.
 */
template <typename Tp>
constexpr Tp inv_mod(Tp a, Tp m) {
  using Sp = std::make_signed_t<Tp>;
  const auto res = euclid_algo<Sp>(a, m).front();
  return static_cast<Tp>(norm_mod<Sp>(res, m));
}
/**
 * \brief Return modular multiplicative inverse of n modulo 2^k where k is the
 * number of digits of Tp.
 *
 * \tparam Tp Unsigned integer required
 * \param n n % 2 == 1 required
 * \return constexpr Tp
 */
template <typename Up>
constexpr Up inv_2p(Up n) {
  static_assert(std::is_unsigned_v<Up>, "Invalid integral type");
  constexpr int Nd = std::numeric_limits<Up>::digits;
  // Let x(k) be the n^(-1) modulo 2^k, then x(2k) = x(k) * (2 - x(k) * n) mod
  // 2^(2k). This is a well known result from newton's method.
  Up ans = 1;  // n^(-1) mod 2^1 = x(2)
  for (int i = 2; i < Nd; i *= 2) { ans *= (2 - ans * n); }
  return ans;
}
/**
 * \brief Find a quadratic non-residue in multiplicative group of Z/nZ. Note
 * that n must be an odd prime number. If the generalized riemann hypothesis is
 * true, the algorithm will run in O(log^2(n))
 *
 * \tparam Tp
 * \param n An 'odd' prime number
 * \return constexpr Tp An integer x which is a quadratic non-residue of
 * multiplicative group modulo n.
 */
template <typename Tp>
constexpr Tp quad_nonresidue(Tp n) {
  // Naive searching. But it won't take much time if GRH is true.
  // If an integer x is quadratic non-residue if and only if x^((n-1)/2)=-1 mod
  // n. Reference: https://en.wikipedia.org/wiki/Euler%27s_criterion

  // n > 2
  Tp ans = 2;
  const Tp exponent = (n - 1) / 2;
  while (pow_mod(ans, exponent, n) != n - 1) { ans += 1; }
  return ans;
}
/**
 * \brief Find one solution for a^x == b (mod m) using Baby-step giant-step
algorithm time complexity: O(sqrt(m))
 *
 * \tparam Tp Integer required
 * \param a 0 <= a < m required
 * \param b 0 <= b < m required
 * \param m m >= 1 required
 * \return std::optional<Tp> An integer x such that a^x == b (mod b).
 */
template <typename Tp>
std::optional<Tp> discrete_log(Tp a, Tp b, Tp m) {
  // sqrt(m) == m, so we need to handle this case.
  if (m == 1) [[unlikely]] { return 0; }
  const Tp n = sqrt(m);
  const Tp an = pow_mod(a, n, m);
  std::map<Tp, Tp> f2;
  for (Tp q = 0, now = b; q <= n; ++q) {
    f2[now] = q;
    now = mul_mod(now, a, m);
  }
  for (Tp p = 1, cur = 1; p <= n; ++p) {
    cur = mul_mod(cur, an, m);
    if (f2.contains(cur)) { return n * p - f2[cur]; }
  }
  return std::nullopt;
}

/**
 * \brief Find a primitive root modulo n or nullopt if it doesn't exist.
primitive root modulo n exists iff:
1. n is 1,2,4 or
2. n is power of an odd prime number (i.e. n=p^k for some prime number p)
or
3. n is of the form 2p^k for some odd prime number p
time complexity: O(Ans*log(phi(n))*log(n)) where phi(x) is euler totient
function
 *
 * \tparam Tp Integer required
 * \param n n >= 1 required
 * \param phi phi = totient(n) required
 * \return constexpr Tp 0 <= res < n will hold.
 */
template <typename Tp>
constexpr auto primitive_root(Tp n, Tp phi) -> std::optional<Tp> {
  assume(n > 0);
  assume(phi > 0);
  if (n == 1) [[unlikely]] { return 0; }
  Tp cur_phi = phi;
  std::array<Tp, 23> factors{};
  int cnt = 0;
  if (cur_phi % 2 == 0) {
    factors[cnt++] = 2;
    do { cur_phi /= 2; } while (cur_phi % 2 == 0);
  }
  if (cur_phi % 3 == 0) {
    factors[cnt++] = 3;
    do { cur_phi /= 3; } while (cur_phi % 3 == 0);
  }
  const Tp sp = sqrt(cur_phi);
  for (Tp i = 5, c = 2; i <= sp; i += c, c ^= 6) {
    if (cur_phi % i == 0) {
      factors[cnt++] = i;
      do { cur_phi /= i; } while (cur_phi % i == 0);
    }
  }
  if (cur_phi > 1) { factors[cnt++] = cur_phi; }
  const auto first = factors.begin();
  const auto last = factors.begin() + cnt;
  for (Tp g = 1; g < n; ++g) {
    if (std::gcd(n, g) != 1) { continue; }
    const auto check = std::find_if(
        first, last, [phi, n, g](Tp p) { return pow_mod(g, phi / p, n) == 1; });
    if (check == last) { return g; }
  }
  return std::nullopt;
}
/**
 * \brief Primitive root modulo prime number n
 *
 * \tparam Tp Integer required
 * \param n prime number required
 * \return constexpr Tp Integer x with 0 <= x < n such that x is a generator of
 * multiplicative group of modulo n.
 */
template <typename Tp>
constexpr Tp primitive_root_prime(Tp n) {
  assume(n >= 1);
  Tp cur_phi = n - 1;
  std::array<Tp, 23> factors{};
  int cnt = 0;
  if (cur_phi % 2 == 0) {
    factors[cnt++] = 2;
    do { cur_phi /= 2; } while (cur_phi % 2 == 0);
  }
  if (cur_phi % 3 == 0) {
    factors[cnt++] = 3;
    do { cur_phi /= 3; } while (cur_phi % 3 == 0);
  }
  const Tp sp = sqrt(cur_phi);
  for (Tp i = 5, c = 2; i <= sp; i += c, c ^= 6) {
    if (cur_phi % i == 0) {
      factors[cnt++] = i;
      do { cur_phi /= i; } while (cur_phi % i == 0);
    }
  }
  if (cur_phi > 1) { factors[cnt++] = cur_phi; }
  const auto first = factors.begin();
  const auto last = factors.begin() + cnt;
  for (Tp g = 1; g < n; ++g) {
    const auto check = std::find_if(
        first, last, [n, g](Tp p) { return pow_mod(g, (n - 1) / p, n) == 1; });
    if (check == last) { return g; }
  }
  detail::unreachable();
}
/**
 * \brief Find one solution for x: x^a == b mod n (where n is prime)
if a == 0 and b == 0, result is undefined
time complexity: O(sqrt(n)+g*log(n)^2) (where g is a primitive root modulo
n)
 *
 * \tparam Tp Integer required
 * \param a
 * \param b 0 <= b < n required.
 * \param n n is prime required
 * \return std::optional<Tp> 0 <= res < n will hold.
 */
template <typename Tp>
std::optional<Tp> discrete_root(Tp a, Tp b, Tp n) {
  const Tp g = primitive_root_prime(n);
  const Tp ga = pow_mod(g, a, n);
  const std::optional<Tp> y = discrete_log(ga, b, n);
  if (y.has_value()) { return pow_mod(g, *y, n); }
  return std::nullopt;
}
}  // namespace algo::detail