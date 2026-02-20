#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
/*
@internal/base/def.hpp
@internal/base/typing.hpp
@internal/base/bit-base.hpp
*/
// makecode
namespace algo {
/**
 * \brief Factorial modulo implementation
 *
 * \tparam MOD
 */
template <detail::modular_integer ModNum>
  requires(ModNum::is_prime_mod())
class factorial {
public:
  /**
   * \brief Precompute factorials up to MOD
   *
   */
  factorial() noexcept = default;
  /**
   * \brief Compute factorials and inverse of them for each integers in [0, n).
   * 0 <= n <= ModNum::mod() must hold.
   *
   * \param n  0 <= n <= MOD should hold
   */
  explicit factorial(int n) noexcept : fac_(n), inv_(n) {
#if !defined(NDEBUG)
    assert(0 <= n && static_cast<uint32_t>(n) <= ModNum::mod());
#endif
    // Edge case
    if (n == 0) [[unlikely]] { return; }
    fac_[0] = 1;
    for (int i = 1; i < n; ++i) { fac_[i] = fac_[i - 1] * ModNum(i); }
    inv_.back() = fac_.back().inv();
    for (int i = n; i-- > 1;) { inv_[i - 1] = inv_[i] * ModNum(i); }
  }
  /**
   * \brief Compute n!.
   *
   * \param n n < limit() must hold.
   * \return ModNum
   */
  ModNum fact(int n) const {
#if !defined(NDEBUG)
    assert(n < limit());
#endif
    return n >= 0 ? fac_[n] : 0;
  }
  /**
   * \brief Return n(n+1)(n+2)...(n+k-1)
   *
   * \param n
   * \param k
   * \return
   */
  ModNum rising_fact(int n, int k) const {
#if !defined(NDEBUG)
    assert(n + k - 1 < limit());
#endif
    return n <= 0 || k < 0 ? 0 : fac_[n + k - 1] * inv_[n - 1];
  }
  /**
   * \brief Return n(n-1)(n-2)...(n-k+1)
   *
   * \param n
   * \param k
   * \return
   */
  ModNum falling_fact(int n, int k) const {
#if !defined(NDEBUG)
    assert(n < limit());
#endif
    return n < 0 || k < 0 || k > n ? 0 : fac_[n] * inv_[n - k];
  }
  // limited inverse factorial
  ModNum inv_fact(int n) const {
#if !defined(NDEBUG)
    assert(n < limit());
#endif
    return n >= 0 ? inv_[n] : 0;
  }
  /**
   * \brief Return binom(n, r)
   *
   * \param n
   * \param r
   * \return mint_type
   */
  ModNum comb(int n, int r) const {
#if !defined(NDEBUG)
    assert(n < limit());
#endif
    return n >= r && r >= 0 ? fac_[n] * inv_[r] * inv_[n - r] : 0;
  }
  /**
   * \brief Return perm(n, r)
   *
   * \param n
   * \param r
   * \return
   */
  ModNum perm(int n, int r) const {
#if !defined(NDEBUG)
    assert(n < limit());
#endif
    return n >= r && r >= 0 ? fac_[n] * inv_[n - r] : 0;
  }
  ModNum catalan(int n) const {
#if !defined(NDEBUG)
    assert(std::max(2 * n, n + 1) < limit());
#endif
    return n < 0 ? 0 : fac_[2 * n] * inv_[n] * inv_[n + 1];
  }
  int limit() const { return static_cast<int>(fac_.size()); }
private:
  std::vector<ModNum> fac_;
  std::vector<ModNum> inv_;
};
}  // namespace algo