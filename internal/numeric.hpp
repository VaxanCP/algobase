#pragma once
#include <array>
// modules
#include "base/bit-base.hpp"
#include "base/def.hpp"
#include "base/math-base.hpp"
#include "base/numeric-base.hpp"

/*
@internal/base/numeric-base.hpp
@internal/base/math-base.hpp
@internal/base/bit-base.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Check if n is prime
 *
 * \tparam Tp Integer required
 * \param n
 * \return true if n is prime
 * \return false otherwise
 */
template <dword_fittable Tp>
constexpr bool prime_test(Tp n) {
  if (n <= 1) { return false; }
  const int s = count_tz(n - 1);
  const Tp d = (n - 1) >> s;
  constexpr std::array<Tp, 4> primes{2, 3, 5, 7};
  for (const Tp p : primes) {
    if (n <= p) { return p == n; }
    // now n > p
    Tp x = pow_mod(p, d, n);
    bool good = x == n - 1 || x == 1;
    for (int r = 1; r < s; ++r) {
      x = mul_mod(x, x, n);
      good |= x == n - 1;
    }
    if (!good) { return false; }
  }
  return true;
}
/**
 * \brief Check if n is prime
 *
 * \tparam Tp Integer required
 * \param n
 * \return true if n is prime
 * \return false otherwise
 */
template <qword_fittable Tp>
constexpr bool prime_test(Tp n) {
  if (n <= 1) { return false; }
  const int s = count_tz(n - 1);
  const Tp d = (n - 1) >> s;
  constexpr std::array<Tp, 12> primes{2,  3,  5,  7,  11, 13,
                                      17, 19, 23, 29, 31, 37};
  for (const Tp p : primes) {
    if (n <= p) { return p == n; }
    // now p < n
    Tp x = pow_mod(p, d, n);
    bool good = x == n - 1 || x == 1;
    for (int r = 1; r < s; ++r) {
      x = mul_mod(x, x, n);
      good |= x == n - 1;
    }
    if (!good) { return false; }
  }
  return true;
}
/**
 * \brief Return value of mobius function at n
 *
 * \tparam Tp Integer required
 * \param n n >= 1 required
 * \return constexpr std::make_signed_t<Tp>
 */
template <typename Tp>
constexpr int mobius(Tp n) {
  assume(n > 0);
  bool parity = false;
  if (n == 1) { return 1; }
  if (n % 2 == 0) {
    if ((n /= 2) % 2 == 0) { return 0; }
    parity ^= true;
  }
  if (n % 3 == 0) {
    if ((n /= 3) % 3 == 0) { return 0; }
    parity ^= true;
  }
  const Tp sn = sqrt(n);
  for (Tp i = 5, c = 2; i <= sn; i += c, c ^= 6) {
    if (n % i == 0) {
      if ((n /= i) % i == 0) { return 0; }
      parity ^= true;
    }
  }
  if (n > 1) { parity ^= 1; }
  return alt(parity);
}
/**
 * \brief Return the value of euler totient function at n
 *
 * \tparam Tp Integer required
 * \param n n >= 1 required
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp totient(Tp n) {
  assume(n > 0);
  Tp ans = n;
  if (n % 2 == 0) {
    ans -= ans / 2;
    do { n /= 2; } while (n % 2 == 0);
  }
  if (n % 3 == 0) {
    ans -= ans / 3;
    do { n /= 3; } while (n % 3 == 0);
  }
  const Tp sn = sqrt(n);
  for (Tp i = 5, c = 2; i <= sn; i += c, c ^= 6) {
    if (n % i == 0) {
      ans -= ans / i;
      do { n /= i; } while (n % i == 0);
    }
  }
  if (n > 1) { ans -= ans / n; }
  return ans;
}
/**
 * \brief Return omega(n)
 *
 * \tparam Tp Integer required
 * \param n n > 0 required
 * \return
 */
template <typename Tp>
constexpr int small_omega(Tp n) {
  assume(n > 0);
  int ans = 0;
  if (n % 2 == 0) {
    ans = ans + 1;
    do { n /= 2; } while (n % 2 == 0);
  }
  if (n % 3 == 0) {
    ans = ans + 1;
    do { n /= 3; } while (n % 3 == 0);
  }
  const Tp sn = sqrt(n);
  for (Tp i = 5, c = 2; i <= sn; i += c, c ^= 6) {
    if (n % i == 0) {
      ans = ans + 1;
      do { n /= i; } while (n % i == 0);
    }
  }
  if (n > 1) { ans = ans + 1; }
  return ans;
}
/**
 * \brief Return Omega(n)
 *
 * \tparam Tp Integer required
 * \param n n > 0 required
 * \return
 */
template <typename Tp>
constexpr int big_omega(Tp n) {
  assume(n > 0);
  int cnt = 0;
  if (n % 2 == 0) {
    do { n /= 2, cnt++; } while (n % 2 == 0);
  }
  if (n % 3 == 0) {
    do { n /= 3, cnt++; } while (n % 3 == 0);
  }
  const Tp sn = sqrt(n);
  for (Tp i = 5, c = 2; i <= sn; i += c, c ^= 6) {
    if (n % i == 0) {
      do { n /= i, cnt++; } while (n % i == 0);
    }
  }
  if (n > 1) { cnt = cnt + 1; }
  return cnt;
}
}  // namespace algo::detail