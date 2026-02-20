#pragma once
#include <algorithm>
#include <cassert>
#include <iterator>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

// lib headers
#include "internal/base/def.hpp"
#include "internal/base/math-base.hpp"
#include "internal/base/typing.hpp"
#include "internal/numeric.hpp"
/*
@internal/numeric.hpp
@internal/base/typing.hpp
@internal/base/math-base.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo {
/**
 * \brief Solve the equation ax+by=gcd(a,b) and return one of the possible
 * solution.
 *
 * \tparam Tp
 * \param a
 * \param b
 * \return Pair of which the first value corresponds to the value of x and
 * the second value corresponds to the gcd(a,b)
 */
template <detail::signed_integer T1, detail::signed_integer T2>
constexpr auto euclid_algo(T1 a, T2 b)
    -> std::array<std::common_type_t<T1, T2>, 3> {
  using Tp = std::common_type_t<T1, T2>;
  return detail::euclid_algo<Tp>(a, b);
}
/**
 * \brief Return the value of euler totient function
 *
 * \tparam Tp integer type
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp totient(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  return detail::totient(n);
}
/**
 * \brief Return the value of mobius function evaluated at n
 *
 * \tparam Tp
 * \param n precondition: n >= 0
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr int mobius(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  return detail::mobius(n);
}
/**
 * \brief Return the number of distinct prime factors of n.
 *
 * \tparam Tp
 * \param n n > 0 required
 * \return
 */
template <detail::integer Tp>
constexpr int small_omega(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  return detail::small_omega(n);
}
/**
 * \brief Return the total numbers of prime factors of n.
 *
 * \tparam Tp
 * \param n precondition: n > 0
 * \return
 */
template <detail::integer Tp>
constexpr int big_omega(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  return detail::big_omega(n);
}
/**
 * \brief Check if given integer is a prime numer or not.
 *
 * \tparam Tp
 * \param n
 * \return true if n is prime
 * \return false otherwise
 */
template <detail::integer Tp>
constexpr bool prime_test(Tp n) {
  return detail::prime_test(n);
}
/**
 * \brief Compute every possible positive quotients of n. i.e Compute [n/i]
 * for 1 <= i <= n. Return array of tuple [l,r,q] denoting [n/i]=q for l <=
 * i <= r.
 *
 * \tparam Tp
 * \param n
 * \return std::vector<std::tuple<Tp, Tp, Tp>>
 */
template <detail::integer Tp>
auto quotients(Tp n) -> std::vector<std::array<Tp, 3>> {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  detail::assume(n > 0);
  std::vector<std::array<Tp, 3>> ans{};
  for (Tp i = 1, j = 1; i <= n; i = j + 1) {
    const Tp q = n / i;
    j = n / q;
    ans.push_back({i, j, q});
  }
  return ans;
}
/**
 * \brief Compute prime factorization of given integer.
 *
 * \tparam Tp
 * \param n n > 0 required.
 * \return Set of prime factors of the given integer and exponents of each. Note
 * that the set is not sorted in increasing order of the factors. The worst cast
 * time complexity is O(sqrt(N)).
 */
template <detail::integer Tp>
auto factorize(Tp n) -> std::vector<std::pair<Tp, int>> {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  detail::assume(n > 0);
  std::vector<std::pair<Tp, int>> ans{};
  ans.reserve(23);
  if (n % 2 == 0) {
    int c = 0;
    do { n /= 2, c++; } while (n % 2 == 0);
    ans.emplace_back(2, c);
  }
  if (n % 3 == 0) {
    int c = 0;
    do { n /= 3, c++; } while (n % 3 == 0);
    ans.emplace_back(3, c);
  }
  const Tp sn = sqrt(n);
  for (Tp i = 5, ec = 2; i <= sn; i += ec, ec ^= 6) {
    if (n % i == 0) {
      int c = 0;
      do { n /= i, c++; } while (n % i == 0);
      ans.emplace_back(i, c);
    }
  }
  if (n > 1) { ans.emplace_back(n, 1); }
  return ans;
}
/**
 * \brief Compute set of divisors of given integer.
 *
 * \tparam Tp
 * \param n n > 0 required.
 * \return Set of divisors of the given integer. Note that the set is not
 * guaranteed to be sorted. Worst cast time complexity is O(sqrt(N)) due to
 * prime factorization.
 */
template <detail::integer Tp>
std::vector<Tp> divisors(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  detail::assume(n > 0);
  std::vector<Tp> res{1};
  for (const auto& [p, c] : factorize(n)) {
    const int sz = static_cast<int>(res.size());
    for (int i = 0; i < sz; ++i) {
      Tp d = res[i];
      for (int j = 0; j < c; ++j) {
        d *= p;
        res.emplace_back(d);
      }
    }
  }
  return res;
}
/**
 * \brief Calculate the values of euler totient function for integers in [0,n).
 * Note that we define phi[0] = 0 for customary. Time complexity is O(NlogN).
 *
 * \tparam Tp
 * \param n n >= 0 required.
 * \return std::vector<int>: The sequence of value of phi
 */
template <detail::integer Tp>
std::vector<Tp> totient_table(Tp n) {
#if !defined(NDEBUG)
  assert(n >= 0);
#endif
  detail::assume(n >= 0);
  std::vector<Tp> phi(n);
  std::iota(phi.begin(), phi.end(), Tp(0));
  // Skip multiples of 2 or 3.
  for (Tp i = 2; i < n; i += 2) { phi[i] -= phi[i] / 2; }
  for (Tp i = 3; i < n; i += 3) { phi[i] -= phi[i] / 3; }
  // This loop go through integer with reminder when divided by 6 is either 1
  // or 5.
  for (Tp i = 5, c = 2; i < n; i += c, c ^= 6) {
    if (phi[i] == i) {
      for (Tp j = i; j < n; j += i) { phi[j] -= phi[j] / i; }
    }
  }
  return phi;
}
/**
 * \brief Calculate the values of mobius function for integers in [0,n). Note
 * that we define mu[0] = 0 for customary. Time complexity is O(NlogN).
 *
 * \tparam Tp
 * \param n n >= 0 required.
 * \return std::vector<int>
 */
template <detail::integer Tp>
std::vector<int> mobius_table(Tp n) {
#if !defined(NDEBUG)
  assert(n >= 0);
#endif
  std::vector<int> mu(n);
  // Edge case
  if (n <= 1) [[unlikely]] { return mu; }
  mu[1] = 1;
  for (Tp i = 2; i < n; ++i) {
    mu[i] -= 1;
    for (Tp j = i + i; j < n; j += i) { mu[j] -= mu[i]; }
  }
  return mu;
}
/**
 * \brief Find all prime numbers in [0,n), and enumerate these integers in
 * increasing order. Time complexity is roughly O(NloglogN).
 *
 * \tparam Tp
 * \param n n >= 0 required
 * \return List of prime numbers less or equal to n.
 */
template <detail::integer Tp>
std::vector<Tp> enumerate_primes(Tp n) {
#if !defined(NDEBUG)
  assert(n >= 0);
#endif
  // Edge case.
  if (n == 0) [[unlikely]] { return {}; }
  const Tp m = (n + 5) / 6;
  // Do sieve for integers in [0, 6*ceil(n/6))
  // 6k+1 <= sqrt(6*m-1) and 6k+5 <= sqrt(6*m-1)
  const Tp kmax = (sqrt(6 * m - 1) + 1) / 6;
  std::vector<bool> is_composite(2 * m);
  is_composite[2 * 0 + 0] = true;
  // This loop goes through all prime numbers up to sqrt(6*m-1)
  for (Tp k = 0; k <= kmax; ++k) {
    if (!is_composite[2 * k + 0]) {
      const Tp pm = 6 * k + 1;
      const Tp s0 = k * (6 * k + 2 * 1) + (1 * 1) / 6;
      const Tp s1 = s0 + 4 * k + 0;
      for (Tp s = s0; s < m; s += pm) { is_composite[2 * s + 0] = true; }
      for (Tp s = s1; s < m; s += pm) { is_composite[2 * s + 1] = true; }
    }
    if (!is_composite[2 * k + 1]) {
      const Tp pm = 6 * k + 5;
      const Tp s0 = k * (6 * k + 2 * 5) + (5 * 5) / 6;
      const Tp s1 = s0 + 2 * k + 1;
      for (Tp s = s0; s < m; s += pm) { is_composite[2 * s + 0] = true; }
      for (Tp s = s1; s < m; s += pm) { is_composite[2 * s + 1] = true; }
    }
  }
  std::vector<Tp> primes{2, 3};
  for (Tp k = 0; k < m; ++k) {
    if (!is_composite[2 * k + 0]) { primes.emplace_back(6 * k + 1); }
    if (!is_composite[2 * k + 1]) { primes.emplace_back(6 * k + 5); }
  }
  while (!primes.empty() && primes.back() >= n) { primes.pop_back(); }
  return primes;
}
}  // namespace algo
