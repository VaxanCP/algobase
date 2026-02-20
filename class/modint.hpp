#pragma once

#include <cassert>
#include <cstdint>
#include <optional>
#include <ostream>
#include <vector>

#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
#include "../internal/modular.hpp"
#include "../internal/numeric.hpp"

/*
@internal/base/typing.hpp
@internal/base/def.hpp
@internal/numeric.hpp
@internal/modular.hpp
*/
// makecode
namespace algo {
/**
 * \brief Modular integer implementation
 *
 * \tparam MOD Modulus
 */
template <uint32_t MOD = 998244353>
  requires(0 < MOD && MOD < (1u << 31))
class modint : detail::modnum_base {
public:
  static constexpr bool IsPrime = detail::prime_test(MOD);
  constexpr modint() noexcept = default;
  template <detail::integer Tp>
  constexpr modint(Tp n) noexcept : value_{transform(n)} {}
  constexpr modint(bool b) noexcept : value_{b ? 1 % MOD : 0} {}
  constexpr modint& operator+=(modint rhs) {
    // Since MOD < 2^31, this->value() + rhs.value() - MOD can safely fit in the
    // signed 32-bit integer.
    value_ += rhs.value() - MOD;
    value_ += (static_cast<int32_t>(value_) < 0) ? MOD : 0;
    return *this;
  }
  constexpr modint& operator-=(modint rhs) {
    value_ -= rhs.value();
    value_ += (static_cast<int32_t>(value_) < 0) ? MOD : 0;
    return *this;
  }
  constexpr modint& operator*=(modint rhs) {
    detail::assume(rhs.value() < MOD);
    detail::assume(value() < MOD);
    value_ = detail::mul_mod(value(), rhs.value(), MOD);
    return *this;
  }
  constexpr modint& operator/=(modint rhs) { return operator*=(rhs.inv()); }
  constexpr modint& operator++() { return operator+=(1); }
  constexpr modint& operator--() { return operator-=(1); }
  constexpr modint operator++(int) {
    const modint res{*this};
    operator++();
    return res;
  }
  constexpr modint operator--(int) {
    const modint res{*this};
    operator--();
    return res;
  }
  /**
   * \brief Unary plus operator.
   *
   * \return constexpr modint
   */
  constexpr modint operator+() const { return *this; }
  /**
   * \brief Unary negation operator.
   *
   * \return constexpr modint
   */
  constexpr modint operator-() const {
    const uint32_t res = value() != 0 ? MOD - value() : 0;
    return from(res);
  }
  constexpr friend auto operator+(modint lhs, modint rhs) -> modint {
    lhs += rhs;
    return lhs;
  }
  constexpr friend auto operator-(modint lhs, modint rhs) -> modint {
    lhs -= rhs;
    return lhs;
  }
  constexpr friend auto operator*(modint lhs, modint rhs) -> modint {
    lhs *= rhs;
    return lhs;
  }
  constexpr friend auto operator/(modint lhs, modint rhs) -> modint {
    lhs /= rhs;
    return lhs;
  }
  constexpr friend bool operator==(modint lhs, modint rhs) {
    return lhs.get() == rhs.get();
  }
  constexpr uint32_t value() const { return value_; }
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
  constexpr modint inv() const
    requires(IsPrime)
  {
#if !defined(NDEBUG)
    assert(value() != 0);
#endif
    return pow(MOD - 2);
  }
  /**
   * \brief Compute multiplicative inverse.
   *
   * \return constexpr modint
   */
  constexpr modint inv() const {
#if !defined(NDEBUG)
    assert(std::gcd(value(), MOD) == 1);
#endif
    // Run extended euclid algorithm.
    return from(detail::inv_mod(value(), MOD));
  }
  /**
   * \brief Compute power to the given exponents.
   *
   * \tparam Tp
   * \param n Exponent. n >= 0 required.
   * \return constexpr modint
   */
  template <detail::integer Tp>
  constexpr modint pow(Tp n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    return detail::binary_exp(*this, n, 1);
  }
  /**
   * \brief Compute modular square root.
   *
   */
  constexpr auto sqrt() const -> std::optional<modint>
    requires(IsPrime)
  {
    if (*this == 0 || *this == 1) { return {*this}; }
    if (pow((MOD - 1) / 2) == -1) { return std::nullopt; }
    constexpr modint Z = detail::quad_nonresidue(MOD);
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
   * \brief
   *
   * \tparam Tp
   * \param x
   * \return constexpr modint
   */
  template <detail::integer Tp>
  constexpr static modint from(Tp x) {
#if !defined(NDEBUG)
    assert(0 <= x && std::cmp_less(x, MOD));
#endif
    return std::bit_cast<modint>(static_cast<uint32_t>(x));
  }
  friend std::ostream& operator<<(std::ostream& os, modint m) {
    os << m.value();
    return os;
  }
  constexpr static uint32_t transform(uint32_t n) { return n % MOD; }
  constexpr static uint32_t transform(int32_t n) {
    n %= static_cast<int32_t>(MOD);
    if (n < 0) { n += static_cast<int32_t>(MOD); }
    return static_cast<uint32_t>(n);
  }
  constexpr static uint32_t transform(uint64_t n) {
    return static_cast<uint32_t>(n % MOD);
  }
  constexpr static uint32_t transform(int64_t n) {
    n %= static_cast<int64_t>(MOD);
    if (n < 0) { n += static_cast<int64_t>(MOD); }
    return static_cast<uint32_t>(n);
  }
  /**
   * \brief Compute modular multiplicative inverse for each integer in [0, n).
   * Note that we define inv[0] = 0 for the customary. Time complexity is O(N).
   *
   * \param n 0 <= n <= MOD must hold.
   * \return inv[i] is the modular inverse of i.
   */
  static std::vector<modint> inv_table(int n)
    requires(IsPrime)
  {
#if !defined(NDEBUG)
    assert(0 <= n && static_cast<uint32_t>(n) <= MOD);
#endif
    std::vector<modint> inv(n);
    // Edge case
    if (n <= 1) [[unlikely]] { return inv; }
    inv[1] = 1;
    for (int i = 2; i < n; ++i) { inv[i] = from(MOD - MOD / i) * inv[MOD % i]; }
    return inv;
  }
private:
  uint32_t value_{0};
};
}  // namespace algo