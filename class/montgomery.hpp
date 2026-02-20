#pragma once

#include <cassert>
#include <cstdint>
#include <optional>
#include <ostream>

#include "../internal/base/bit-base.hpp"
#include "../internal/base/numeric-base.hpp"
#include "../internal/base/typing.hpp"
#include "../internal/modular.hpp"
#include "../internal/numeric.hpp"

/*
@internal/modular.hpp
@internal/base/typing.hpp
@internal/base/bit-base.hpp
@internal/base/numeric-base.hpp
@internal/numeric.hpp
*/

// makecode
namespace algo {
/**
 * \brief Montgomery modint implementation. Montgomery multiplication is an
 * optimization technique involving modular multiplication which is widely used
 * in cryptography. For integers x and y, x is said to be transformed into y in
 * montgomery space if and only if:
 *  1. y mod N = x*2^32 mod N
 *  2. 0 <= y < MOD
 *
 * \tparam MOD
 */
template <uint32_t MOD = 998244353>
  requires(MOD % 2 == 1 && MOD < (1u << 31))
class montgomery_modint : detail::modnum_base {
public:
  static constexpr uint32_t R1 = (-MOD) % MOD;                       // 2^32
  static constexpr uint32_t R2 = -static_cast<uint64_t>(MOD) % MOD;  // 2^64
  static constexpr uint32_t R3 = detail::mul_mod(R1, R2, MOD);       // 2^96
  static constexpr uint32_t Nd = detail::inv_2p(MOD);                // inv(MOD)
  // Check if MOD is a prime number or not (Which is needed to speed up the
  // inverse computation)
  static constexpr bool IsPrime = detail::prime_test(MOD);
  constexpr montgomery_modint() noexcept = default;
  template <detail::integer Tp>
  constexpr montgomery_modint(Tp n) noexcept : value_{transform(n)} {}
  constexpr montgomery_modint(bool b) noexcept : value_{b ? R1 : 0} {}
  constexpr auto operator+=(montgomery_modint rhs) -> montgomery_modint& {
    // Keep this->value_ inside the range [0, M)
    value_ += rhs.get() - MOD;
    // If value_ is negative, add M
    value_ += (static_cast<int32_t>(value_) >> 31) & MOD;
    return *this;
  }
  constexpr auto operator-=(montgomery_modint rhs) -> montgomery_modint& {
    // Keep this->value_ inside the range [0, M)
    value_ -= rhs.get();
    // If value_ is negative, add M
    value_ += (static_cast<int32_t>(value_) >> 31) & MOD;
    return *this;
  }
  constexpr auto operator*=(montgomery_modint rhs) -> montgomery_modint& {
    value_ = mult(get(), rhs.get());
    return *this;
  }
  constexpr auto operator/=(montgomery_modint rhs) -> montgomery_modint& {
    // a/b = a*b^(-1)
    return operator*=(rhs.inv());
  }
  constexpr auto operator++() -> montgomery_modint& {
    return operator+=(from(R1));
  }
  constexpr auto operator--() -> montgomery_modint& {
    return operator-=(from(R1));
  }
  constexpr auto operator++(int) -> montgomery_modint {
    const montgomery_modint res{*this};
    operator++();
    return res;
  }
  constexpr auto operator--(int) -> montgomery_modint {
    const montgomery_modint res{*this};
    operator--();
    return res;
  }
  constexpr auto operator+() const -> montgomery_modint { return *this; }
  constexpr auto operator-() const -> montgomery_modint {
    // To keep the value inside [0,M)
    const uint32_t res = get() != 0 ? MOD - get() : 0;
    return from(res);
  }
  constexpr friend auto operator+(montgomery_modint lhs, montgomery_modint rhs)
      -> montgomery_modint {
    lhs += rhs;
    return lhs;
  }
  constexpr friend auto operator-(montgomery_modint lhs, montgomery_modint rhs)
      -> montgomery_modint {
    lhs -= rhs;
    return lhs;
  }
  constexpr friend auto operator*(montgomery_modint lhs, montgomery_modint rhs)
      -> montgomery_modint {
    lhs *= rhs;
    return lhs;
  }
  constexpr friend auto operator/(montgomery_modint lhs, montgomery_modint rhs)
      -> montgomery_modint {
    lhs /= rhs;
    return lhs;
  }
  constexpr friend bool operator==(montgomery_modint lhs,
                                   montgomery_modint rhs) {
    return lhs.get() == rhs.get();
  }
  /**
   * \brief Get reduced value.
   *
   * \return
   */
  constexpr uint32_t value() const { return reduce(value_); }
  /**
   * \brief Get raw value
   *
   * \return
   */
  constexpr uint32_t get() const { return value_; }
  static constexpr uint32_t mod() { return MOD; }
  /**
   * \brief Check if the modulus is a prime number or not.
   *
   * \return true
   * \return false
   */
  static constexpr bool is_prime_mod() { return IsPrime; }
  /**
   * \brief Compute multiplicative inverse.
   *
   */
  constexpr auto inv() const -> montgomery_modint
    requires(IsPrime)
  {
#if !defined(NDEBUG)
    assert(get() != 0);
#endif
    return pow(MOD - 2);
  }
  /**
   * \brief Compute modular multiplicative inverse.
   *
   * \return montgomery_modint
   */
  constexpr auto inv() const -> montgomery_modint {
    // This will be (Rx)^-1, so we will multiply it by R^3.
    const uint32_t res = detail::inv_mod(get(), MOD);
    return from(mult(res, R3));
  }
  /**
   * \brief Compute modular exponentiation.
   *
   * \tparam Tp
   * \param n
   * \return montgomery_modint
   */
  template <detail::integer Tp>
  constexpr auto pow(Tp n) const -> montgomery_modint {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    return detail::binary_exp(*this, n, 1);
  }
  /**
   * \brief Compute square root of this value.
   *
   */
  constexpr auto sqrt() const -> std::optional<montgomery_modint>
    requires(IsPrime)
  {
    if (*this == 0 || *this == 1) { return {*this}; }
    if (pow((MOD - 1) / 2) == -1) { return std::nullopt; }
    constexpr montgomery_modint Z = detail::quad_nonresidue(MOD);
    constexpr uint32_t K = detail::count_tz(MOD - 1);
    constexpr uint32_t C = (MOD - 1) >> K;
    auto res = pow((C + 1) / 2);
    auto cur = pow(C);
    auto fac = Z.pow(C);
    for (uint32_t i = K - 1; i > 0 && cur != 1; --i) {
      const auto nfac = fac * fac;
      if (cur.pow(1u << (i - 1)) == -1) {
        res *= fac;
        cur *= nfac;
      }
      fac = nfac;
    }
    return {res};
  }
  /**
   * \brief Given an integer n which is already in the montgomery space, cast it
   * to the montgomery_modint. Note that this is bit-wise cast and doesn't incur
   * any latency.
   *
   * \tparam Tp
   * \param n
   * \return constexpr montgomery_modint
   */
  template <detail::integer Tp>
  static constexpr auto from(Tp n) -> montgomery_modint {
#if !defined(NDEBUG)
    assert(0 <= n && std::cmp_less(n, MOD));
#endif
    return std::bit_cast<montgomery_modint>(static_cast<uint32_t>(n));
  }
  friend std::ostream& operator<<(std::ostream& os, montgomery_modint m) {
    os << m.value();
    return os;
  }
  // Given an integer which is in montgomery space, return the original residue
  // class.
  static constexpr uint32_t reduce(uint32_t n) {
    const uint32_t q = n * Nd;
    const uint64_t r = detail::safe_mul(q, MOD);
    const uint32_t g = static_cast<uint32_t>(r >> 32);
    const uint32_t d = 0 - g;
    // d may be negative
    return d + ((static_cast<int32_t>(d) >> 31) & MOD);
  }
  // Given two integers in montgomery space, multiply them.
  // The result is also in montgomery space.
  static constexpr uint32_t mult(uint32_t x, uint32_t y) {
    const uint64_t prod = detail::safe_mul(x, y);
    const uint32_t l = static_cast<uint32_t>(prod);
    const uint32_t h = static_cast<uint32_t>(prod >> 32);
    const uint32_t q = l * Nd;
    const uint64_t r = detail::safe_mul(q, MOD);
    const uint32_t g = static_cast<uint32_t>(r >> 32);
    const uint32_t d = h - g;
    return d + ((static_cast<int32_t>(d) >> 31) & MOD);
  }
  // Given an integer, transform it into the montgomery space.
  static constexpr uint32_t transform(uint32_t n) { return mult(n, R2); }
  // Given an integer, transform it into the montgomery space.
  static constexpr uint32_t transform(int32_t n) {
    n %= static_cast<int32_t>(MOD);
    if (n < 0) { n += static_cast<int32_t>(MOD); }
    return mult(static_cast<uint32_t>(n), R2);
  }
  // Given an integer, transform it into the montgomery space.
  static constexpr uint32_t transform(uint64_t n) {
    return mult(static_cast<uint32_t>(n % MOD), R2);
  }
  // Given an integer, transform it into the montgomery space.
  static constexpr uint32_t transform(int64_t n) {
    n %= static_cast<int64_t>(MOD);
    if (n < 0) { n += static_cast<int64_t>(MOD); }
    return mult(static_cast<uint32_t>(n), R2);
  }
  /**
   * \brief Compute modular multiplicative inverse for each integers in [0, n).
   * Note that we define inv[0] = 0 for the customary.
   * \param n 0 <= n <= MOD must hold.
   *
   */
  static auto inv_table(int n) -> std::vector<montgomery_modint>
    requires(IsPrime)
  {
#if !defined(NDEBUG)
    assert(0 <= n && static_cast<uint32_t>(n) <= MOD);
#endif
    std::vector<montgomery_modint> inv(n);
    // Edge case
    if (n <= 1) [[unlikely]] { return inv; }
    inv[1] = 1;
    for (int i = 2; i < n; ++i) {
      inv[i] = montgomery_modint(MOD - MOD / i) * inv[MOD % i];
    }
    return inv;
  }
private:
  uint32_t value_{0};
};
}  // namespace algo