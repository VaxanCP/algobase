#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/fft.hpp"
#include "./modint.hpp"
#include "./montgomery.hpp"

/*
@internal/base/def.hpp
@internal/base/bit-base.hpp
@internal/fft.hpp
@class/modint.hpp
@class/montgomery.hpp
*/
// makecode
namespace algo {
namespace detail {
template <uint32_t MOD>
std::vector<modint<MOD>> convolution1(const std::vector<modint<MOD>>& f1,
                                      const std::vector<modint<MOD>>& f2,
                                      int n) noexcept {
  const size_t n1 = f1.size();
  const size_t n2 = f2.size();
  const size_t block_sz = size_t(n) << 1;
  const uint32_t* c1 = reinterpret_cast<const uint32_t*>(f1.data());
  const uint32_t* c2 = reinterpret_cast<const uint32_t*>(f2.data());
  uint32_t* buffer = fft::get_storage(block_sz);
  uint32_t* p1 = buffer + 0 * n;
  uint32_t* p2 = buffer + 1 * n;
  memcpy(p1, c1, n1 * sizeof(uint32_t));
  memcpy(p2, c2, n2 * sizeof(uint32_t));
  fft::convolution<MOD>(p1, p2, n);
  modint<MOD>* out = reinterpret_cast<modint<MOD>*>(p1);
  std::vector<modint<MOD>> res(out, out + n);
  fft::free_storage(buffer);
  return res;
}
template <uint32_t MOD>
std::vector<modint<MOD>> convolution3(const std::vector<modint<MOD>>& f1,
                                      const std::vector<modint<MOD>>& f2,
                                      int n) noexcept {
  constexpr uint32_t M1 = 167772161;
  constexpr uint32_t M2 = 469762049;
  constexpr uint32_t M3 = 754974721;
  const size_t n1 = f1.size();
  const size_t n2 = f2.size();
  const size_t block_sz = size_t(n) << 2;
  const uint32_t* c1 = reinterpret_cast<const uint32_t*>(f1.data());
  const uint32_t* c2 = reinterpret_cast<const uint32_t*>(f2.data());
  uint32_t* buffer = fft::get_storage(block_sz);
  uint32_t* p1 = buffer + 0 * n;
  uint32_t* p2 = buffer + 1 * n;
  uint32_t* p3 = buffer + 2 * n;
  uint32_t* p4 = buffer + 3 * n;
  // first
  memcpy(p1, c1, n1 * sizeof(uint32_t));
  memcpy(p4, c2, n2 * sizeof(uint32_t));
  fft::convolution<M1>(p1, p4, n);
  // second
  memcpy(p2, c1, n1 * sizeof(uint32_t));
  memcpy(p4, c2, n2 * sizeof(uint32_t));
  fft::convolution<M2>(p2, p4, n);
  // third
  memcpy(p3, c1, n1 * sizeof(uint32_t));
  memcpy(p4, c2, n2 * sizeof(uint32_t));
  fft::convolution<M3>(p3, p4, n);
  // final
  std::vector<modint<MOD>> res(n);
  uint32_t* out = reinterpret_cast<uint32_t*>(res.data());
  fft::consolidate<MOD>(p1, p2, p3, out, n);
  fft::free_storage(buffer);
  return res;
}
}  // namespace detail
/**
 * \brief Formal power series implementation
 *
 * \tparam MOD
 */
template <uint32_t MOD = 998244353>
class poly {
  static constexpr bool is_large(int n) { return n >= 128; }
public:
  using mint_type = modint<MOD>;
  poly() = default;
  template <std::input_iterator IIter>
    requires std::convertible_to<std::iter_value_t<IIter>, mint_type>
  poly(IIter first, IIter last) noexcept : c_(first, last) {
    drop_zeros();
  }
  poly(std::vector<mint_type> vec) noexcept : c_(std::move(vec)) {
    drop_zeros();
  }
  template <detail::integer Tp>
  poly(std::initializer_list<Tp> lst) noexcept : c_(lst.begin(), lst.end()) {
    drop_zeros();
  }
  poly(std::initializer_list<mint_type> lst) noexcept : c_(lst) {
    drop_zeros();
  }
  poly& operator+=(const poly& rhs) {
    if (size() < rhs.size()) { resize(rhs.c_.size()); }
    for (int i = 0; i < rhs.size(); ++i) { c_[i] += rhs.c_[i]; }
    drop_zeros();
    return *this;
  }
  poly& operator+=(mint_type x) {
    if (is_zero()) { c_.emplace_back(0); }
    c_[0] += x;
    drop_zeros();
    return *this;
  }
  poly& operator-=(const poly& rhs) {
    if (size() < rhs.size()) { resize(rhs.c_.size()); }
    for (int i = 0; i < rhs.size(); ++i) { c_[i] -= rhs.c_[i]; }
    drop_zeros();
    return *this;
  }
  poly& operator-=(mint_type x) {
    if (is_zero()) { c_.emplace_back(0); }
    c_[0] -= x;
    drop_zeros();
    return *this;
  }
  poly& operator*=(const poly& rhs) {
    *this = (*this) * rhs;
    return *this;
  }
  poly& operator*=(mint_type x) {
    for (auto& c : c_) { c *= x; }
    drop_zeros();
    return *this;
  }
  poly& operator/=(const poly& rhs) {
    *this = *this / rhs;
    return *this;
  }
  poly& operator%=(const poly& rhs) {
    *this = *this % rhs;
    return *this;
  }
  poly operator*(const poly& rhs) const {
    const int m = std::min(size(), rhs.size());
    return is_large(m) ? fast_mul(rhs) : slow_mul(rhs);
  }
  poly operator/(const poly& rhs) const {
#if !defined(NDEBUG)
    assert(!rhs.is_zero());
#endif
    const int m = deg() - rhs.deg() + 1;
    return is_large(m) ? fast_div(rhs) : slow_div(rhs);
  }
  poly operator%(const poly& rhs) const {
#if !defined(NDEBUG)
    assert(!rhs.is_zero());
#endif
    const poly D = *this / rhs;
    return *this - D * rhs;
  }
  friend poly operator-(poly p) {
    for (auto& x : p.c_) { x = -x; }
    return p;
  }
  friend poly operator+(poly lhs, const poly& rhs) {
    lhs += rhs;
    return lhs;
  }
  friend poly operator-(poly lhs, const poly& rhs) {
    lhs -= rhs;
    return lhs;
  }
  // this one enables something like: poly + int
  friend poly operator+(poly lhs, mint_type rhs) {
    lhs += rhs;
    return lhs;
  }
  friend poly operator+(mint_type lhs, poly rhs) {
    rhs += lhs;
    return rhs;
  }
  friend poly operator-(poly lhs, mint_type rhs) {
    lhs -= rhs;
    return lhs;
  }
  friend poly operator-(mint_type lhs, poly rhs) {
    rhs -= lhs;
    return -rhs;
  }
  friend poly operator*(poly lhs, mint_type rhs) {
    lhs *= rhs;
    return lhs;
  }
  friend poly operator*(mint_type lhs, poly rhs) {
    rhs *= lhs;
    return rhs;
  }
  mint_type operator[](int i) const {
    if (i < 0 || i > deg()) { return 0; }
    return c_[i];
  }
  mint_type get_unsafe(int i) const { return c_[i]; }
  mint_type operator()(mint_type x) const {
    mint_type ans = 0;
    for (const auto c : c_) {
      ans *= x;
      ans += c;
    }
    return ans;
  }
  /**
   * \brief Take only first n terms
   *
   * \param n
   * \return poly
   */
  poly take(int n) const& {
#if !defined(NDEBUG)
    assert(n > 0);
#endif
    return std::vector(c_.begin(), c_.begin() + std::min(n, size()));
  }
  poly take(int n) && {
#if !defined(NDEBUG)
    assert(n > 0);
#endif
    resize(std::min(n, size()));
    return std::move(c_);
  }
  /**
   * \brief Return p(x) * x^n
   *
   * \param n
   * \return poly
   */
  poly shift(int n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    std::vector<mint_type> out(n + size());
    std::copy_backward(c_.begin(), c_.end(), out.end());
    return out;
  }
  /**
   * \brief Return f(x) / x^n
   *
   * \param n
   * \return poly
   */
  poly drop(int n) const& {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    return std::vector(c_.begin() + std::min(n, size()), c_.end());
  }
  poly drop(int n) && {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    std::shift_left(c_.begin(), c_.end(), n);
    resize(std::max(0, size() - n));
    return std::move(c_);
  }
  /**
   * \brief Return f(1/x) * x^n
   *
   * \param n
   * \return poly
   */
  poly rev(int n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    std::vector<mint_type> res(std::max(size(), n + 1));
    std::copy_backward(c_.rbegin(), c_.rend(), res.end());
    return res;
  }
  /**
   * \brief Compute derivative of p(x)
   *
   * \return poly
   */
  poly diff() const& {
    std::vector<mint_type> res(std::max(0, deg()));
    for (int i = 1; i < size(); ++i) { res[i - 1] = c_[i] * i; }
    return res;
  }
  poly diff() && {
    const int nsz = std::max(0, size() - 1);
    for (int i = 0; i < nsz; ++i) { c_[i] = c_[i + 1] * mint_type(i + 1); }
    resize(nsz);
    return std::move(c_);
  }
  /**
   * \brief Compute integral(p(x))
   *
   * \return poly
   */
  poly integrate() const& {
    std::vector<mint_type> res(size() + 1);
    for (int i = 0; i <= deg(); ++i) { res[i + 1] = c_[i] / mint_type(i + 1); }
    return res;
  }
  poly integrate() && {
    c_.emplace_back(0);
    for (int i = size() - 1; i > 0; --i) { c_[i] = c_[i - 1] / mint_type(i); }
    c_[0] = 0;
    return std::move(c_);
  }

  /**
   * \brief Return F(x)^-1 mod x^n
   *
   * \param n
   * \return poly
   */
  poly inv(int n) const {
#if !defined(NDEBUG)
    assert(size() && c_[0] != 0);
    assert(n > 0);
#endif
    int pw = 1;
    poly ans = {c_[0].inv()};
    while (pw < n) {
      poly C = ans * take(pw << 1);
      poly D = ans * std::move(C).take(pw << 1).drop(pw);
      ans -= std::move(D).take(pw).shift(pw);
      pw <<= 1;
    }
    return std::move(ans).take(n);
  }
  /**
   * \brief Return p(x)^k mod x^n in O(n logn long k)
   *
   * \param k
   * \param n
   * \return poly
   */
  poly binpow(int k, int n) const {
#if !defined(NDEBUG)
    assert(k >= 0);
    assert(n > 0);
#endif
    poly res = {mint_type::from(1)};
    poly mul = take(n);
    while (k > 0) {
      if (k & 0x01) { res = (res * mul).take(n); }
      mul = (mul * mul).take(n);
      k >>= 1;
    }
    return res;
  }
  /**
   * \brief Compute log(p(x)) mod x^n
   *
   * \param n
   * \return poly
   */
  poly log(int n) const {
#if !defined(NDEBUG)
    assert(n > 0);
    assert(size() && c_[0] != 0);
#endif
    return (diff().take(n) * inv(n)).integrate().take(n);
  }
  /**
   * \brief Compute exp(p(x)) mod x^n
   *
   * \param n
   * \return poly
   */
  poly exp(int n) const {
#if !defined(NDEBUG)
    assert(n > 0);
    assert(size() == 0 || c_[0] == 0);
#endif
    int pw = 1;
    poly ans = {mint_type::from(1)};
    while (pw < n) {
      const poly H = take(pw << 1).drop(pw) - ans.log(pw << 1).drop(pw);
      ans += (ans * H).take(pw).shift(pw);
      pw <<= 1;
    }
    return std::move(ans).take(n);
  }
  mint_type lead() const { return is_zero() ? 0 : c_.back(); }
  std::vector<mint_type> data() const& { return c_; }
  std::vector<mint_type> data() && { return std::move(c_); }
  int deg() const { return size() - 1; }
  int size() const { return static_cast<int>(c_.size()); }
  bool is_zero() const { return c_.empty(); }
private:
  void drop_zeros() noexcept {
    while (!c_.empty() && c_.back() == 0) { c_.pop_back(); }
  }
  void resize(size_t n) noexcept {
    detail::assume(c_.size() < c_.max_size());
    c_.resize(n);
  }
  poly slow_mul(const poly& rhs) const {
    if (!is_zero() && !rhs.is_zero()) [[likely]] {
      const int n1 = size();
      const int n2 = rhs.size();
      const int n = n1 + n2 - 1;
      std::vector<mint_type> res(n);
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) { res[i + j] += c_[i] * rhs.c_[j]; }
      }
      return res;
    } else {
      return {};
    }
  }
  poly fast_mul(const poly& rhs) const {
    constexpr int Cap = 1 << detail::count_tz(MOD - 1);
    constexpr int Lim = 1 << 23;
    const int n1 = size();
    const int n2 = rhs.size();
    const int n = detail::ceil_pow2(n1 + n2 - 1);
    if (Cap >= Lim || n <= Cap) {
      return detail::convolution1(c_, rhs.c_, n);
    } else {
      return detail::convolution3(c_, rhs.c_, n);
    }
  }
  poly slow_div(const poly& rhs) const {
    const int d1 = deg();
    const int d2 = rhs.deg();
    const int m = d1 - d2 + 1;
    if (d1 < d2) { return poly{}; }
    std::vector<mint_type> out(m), temp(m);
    std::copy(c_.begin() + d2, c_.end(), temp.begin());
    const mint_type x = rhs.c_.back().inv();
    for (int i = m - 1; i >= 0; --i) {
      out[i] = temp[i];
      if (temp[i] != 0) {
        out[i] *= x;
        const int q = std::max(d2 - i, 0);
        for (int j = d2; j >= q; --j) {
          temp[i + j - d2] -= out[i] * rhs.c_[j];
        }
      }
    }
    return out;
  }
  // precondition: deg() >>> rhs.deg()
  poly fast_div(const poly& rhs) const {
    const int d1 = deg();
    const int d2 = rhs.deg();
    const int m = d1 - d2 + 1;
    return (rev(d1).take(m) * rhs.rev(d2).inv(m)).take(m).rev(m - 1);
  }
  std::vector<mint_type> c_;
};
}  // namespace algo