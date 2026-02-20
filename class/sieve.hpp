#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "../internal/base/def.hpp"
#include "../internal/base/math-base.hpp"

/*
@internal/base/math-base.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo {
class sieve {
public:
  /**
   * \brief Default constructor. Create empty sieve engine.
   *
   */
  sieve() noexcept = default;
  /**
   * \brief Construct a new sieve engine. Precompute a prime factor for each
   * integers in [0, n). Note that we don't define prime factors for 0.
   *
   * \param n A integer denoting the upper bound.
   */
  explicit sieve(int n) noexcept : lim_{n}, factors_(2 * ((n + 5) / 6)) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    // Edge case
    if (n == 0) [[unlikely]] { return; }
    /*
      Do sieve for integers in [0, 6*ceil(n/6))
    */
    // 6 * k + 1 => 2 * k + 0
    // 6 * k + 5 => 2 * k + 1
    // f[i] = i / 3 mapped
    const int m = (n + 5) / 6;
    const int kmax = (detail::sqrt(6 * m - 1) + 1) / 6;
    // Let m = ceil(n / 6), then this loop covers every primes up to 6*m-1.
    // If n is multiple of 6, then it can be ignored because n can never be a
    // prime number.
    // factor[1] = 1
    factors_[2 * 0 + 0] = 1;
    for (int k = 0; k <= kmax; ++k) {
      if (factors_[2 * k + 0] == 0) {
        const int pm = 6 * k + 1;
        const int s0 = k * (6 * k + 2 * 1) + (1 * 1) / 6;
        const int s1 = s0 + 4 * k + 0;
        for (int s = s0; s < m; s += pm) { factors_[2 * s + 0] = pm; }
        for (int s = s1; s < m; s += pm) { factors_[2 * s + 1] = pm; }
      }
      if (factors_[2 * k + 1] == 0) {
        const int pm = 6 * k + 5;
        const int s0 = k * (6 * k + 2 * 5) + (5 * 5) / 6;
        const int s1 = s0 + 2 * k + 1;
        for (int s = s0; s < m; s += pm) { factors_[2 * s + 0] = pm; }
        for (int s = s1; s < m; s += pm) { factors_[2 * s + 1] = pm; }
      }
    }
  }
  /**
   * \brief Compute the prime factorization of given integer. Note that the
   * factors are not sorted in increasing order.
   *
   * \param n The integer to factorize. 0 < n < limit() must hold.
   * \return std::vector<std::array<int, 2>> Prime factorization of given
   * integer. For each element p, p[0] denotes the prime factor and p[1] denotes
   * the exponents of it.
   */
  auto factorize(int n) const -> std::vector<std::array<int, 2>> {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    std::vector<std::array<int, 2>> ans;
    ans.reserve(20);
    while (n > 1) {
      const int p = get_factor(n);
      int c = 0;
      do { n /= p, c++; } while (n % p == 0);
      ans.push_back({p, c});
    }
    return ans;
  }
  /**
   * \brief Compute the set of divisors of a given integer. Note that the
   * divisors are not guaranteed to be sorted and it includes the trivial
   * divisors (i.e., 1 and the given integer).
   *
   * \param n The integer of which the set of divisors are computed. 0 < n <
   * limit() must hold. \return std::vector<int> Set of divisors of given
   * integer.
   */
  std::vector<int> divisors(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    std::vector<int> ans{1};
    while (n > 1) {
      const int p = get_factor(n);
      int c = 0;
      do { n /= p, c++; } while (n % p == 0);
      const int sz = static_cast<int>(ans.size());
      for (int i = 0; i < sz; ++i) {
        for (int j = 0, d = ans[i]; j < c; ++j) {
          d *= p;
          ans.emplace_back(d);
        }
      }
    }
    return ans;
  }
  /**
   * \brief Compute the value of euler totient function at given integer.
   * \param n A integer. 0 < n < limit() should hold.
   * \return int The value of euler totient function at given integer.
   */
  int totient(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    int ans = n;
    while (n > 1) {
      const int p = get_factor(n);
      ans -= ans / p;
      do { n /= p; } while (n % p == 0);
    }
    return ans;
  }
  /**
   * \brief Compute the value of mobius function at given integer.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int The value of mobius function at given integer.
   */
  int mobius(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    bool parity = false;
    while (n > 1) {
      const int p = get_factor(n);
      if ((n /= p) % p == 0) { return 0; }
      parity ^= 1;
    }
    return detail::alt(parity);
  }
  /**
   * \brief Compute the number of distinct prime factors of given integer.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int The number of distinct prime factors of given integer.
   */
  int small_omega(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    int ans = 0;
    while (n > 1) {
      const int p = get_factor(n);
      ans = ans + 1;
      do { n /= p; } while (n % p == 0);
    }
    return ans;
  }
  /**
   * \brief Compute the number of total prime factors of given integer.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int The number of total prime factors of given integer.
   */
  int big_omega(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    int ans = 0;
    while (n > 1) {
      const int p = get_factor(n);
      do { n /= p, ans++; } while (n % p == 0);
    }
    return ans;
  }
  /**
   * \brief Compute the radical number of given integer.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int Radical number of given integer.
   */
  int radical(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    int ans = 1;
    while (n > 1) {
      const int p = get_factor(n);
      ans *= p;
      do { n /= p; } while (n % p == 0);
    }
    return ans;
  }
  /**
   * \brief Compute the number of divisors of a given integer.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int The number of divisors, which includes 1 and the given integer
   * itself.
   */
  int sigma(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    int ans = 1;
    while (n > 1) {
      const int p = get_factor(n);
      int c = 0;
      do { n /= p, c++; } while (n % p == 0);
      ans *= 1 + c;
    }
    return ans;
  }
  /**
   * \brief Check if a given integer is a prime number.
   *
   * \param n A integer to check.
   * \return true if the given integer is a prime number, false otherwise.
   */
  bool is_prime(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    return n == get_factor(n);
  }
  /**
   * \brief Given an integer, compute one of its prime factor (not necessalilly
   * the smallest one). If the integer is a prime number, return as is.
   *
   * \param n A integer. 0 < n < limit() should hold.
   * \return int A prime factor of given integer.
   */
  int get_factor(int n) const {
#if !defined(NDEBUG)
    assert(0 < n && n < limit());
#endif
    detail::assume(0 < n);
    /// Trivial cases
    if (n % 2 == 0) { return 2; }
    if (n % 3 == 0) { return 3; }
    return get_nontrivial_factor(n);
  }
  /**
   * \brief Upper bound
   *
   * \return
   */
  int limit() const { return lim_; }
private:
  int get_nontrivial_factor(int n) const {
    detail::assume(0 < n);
    // n = 6k+1 => n/3 = 2k+0
    // n = 6k+5 => n/3 = 2k+1
    const int index = n / 3;
    return factors_[index] ? factors_[index] : n;
  }
  int lim_;
  std::vector<int> factors_;
};
}  // namespace algo