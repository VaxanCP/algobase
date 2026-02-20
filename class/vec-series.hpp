#pragma once

#include <immintrin.h>

#include <algorithm>
#include <cassert>
#include <span>
#include <utility>
#include <vector>

#include "../internal/numeric.hpp"
#include "../internal/simd/fft.hpp"
#include "../internal/simd/helpers.hpp"
#include "../internal/simd/modint.hpp"
#include "./montgomery.hpp"

/*
@internal/simd/fft.hpp
@internal/simd/helpers.hpp
@internal/simd/modint.hpp
@internal/numeric.hpp
@class/montgomery.hpp
*/

// makecode
namespace algo {
/**
 * \brief Super fast FPS implementation dedicated for NTT friendly modulus. MOD
 * must be a prime number ans MOD % 8 == 1 must hold. Like in fast_fps
 * implementation, we assume that at any time size() will not exceed
 * 2^(count_tz(MOD-1)). The internal representation of coefficients models
 * montgomery_modint<MOD>.
 *
 * \tparam MOD
 */
template <uint32_t MOD = 998244353>
  requires(detail::prime_test(MOD) && MOD % 8 == 1)
class vec_series {
  using mint_type = montgomery_modint<MOD>;
  using mintx8_type = detail::simd::montgomery_modintx8<MOD>;
public:
  /**
   * \brief Default constructor.
   *
   */
  vec_series() noexcept = default;
  /**
   * \brief Construct new vec_series with specified coefficients.
   *
   * \param vec
   */
  vec_series(const std::vector<mint_type>& vec) noexcept
      : c_((vec.size() + 7) / 8) {
    const int n = static_cast<int>(vec.size());
    detail::assume(n >= 0);
    const int n8 = n / 8;
    // Temporary buffer
    alignas(__m256i) std::array<mint_type, 8> buf;
    for (int i = 0; i < n8; ++i) {
      for (int j = 0; j < 8; ++j) { buf[j] = vec[i * 8 + j]; }
      c_[i].load_aligned(buf.data());
    }
    if (n % 8 != 0) {
      for (int j = 0; j < n % 8; ++j) { buf[j] = vec[n8 * 8 + j]; }
      for (int j = n % 8; j < 8; ++j) { buf[j] = 0; }
      c_[n8].load_aligned(buf.data());
    }
  }
  /**
   * \brief Construct new vec_series from the given elements.
   *
   * \param lst
   */
  vec_series(std::initializer_list<mint_type> lst) noexcept
      : c_((lst.size() + 7) / 8) {
    const int n = static_cast<int>(lst.size());
    detail::assume(n >= 0);
    const int n8 = n / 8;
    // Temporary buffer
    alignas(__m256i) std::array<mint_type, 8> buf;
    for (int i = 0; i < n8; ++i) {
      std::copy_n(lst.begin() + 8 * i, 8, buf.begin());
      c_[i].load_aligned(buf.data());
    }
    if (n % 8 != 0) {
      buf.fill(0);
      std::copy_n(lst.begin() + 8 * n8, n % 8, buf.begin());
      c_[n8].load_aligned(buf.data());
    }
  }
  /**
   * \brief Construct new vec_series
   *
   * \param vec
   */
  vec_series(std::vector<mintx8_type> vec) noexcept : c_{std::move(vec)} {}
  /**
   * \brief Add rhs to this series.
   *
   * \param rhs
   * \return vec_series&
   */
  vec_series& operator+=(const vec_series& rhs) {
    // Resize to the larger one
    if (c_.size() < rhs.c_.size()) { c_.resize(rhs.c_.size()); }
    // Do addition
    for (int i = 0; i < rhs.size_8(); ++i) { c_[i] += rhs.c_[i]; }
    return *this;
  }
  /**
   * \brief Subtract rhs from this series.
   *
   * \param rhs
   * \return vec_series&
   */
  vec_series& operator-=(const vec_series& rhs) {
    // Resize to the larger one
    if (c_.size() < rhs.c_.size()) { c_.resize(rhs.c_.size()); }
    // Do subtraction
    for (int i = 0; i < rhs.size_8(); ++i) { c_[i] -= rhs.c_[i]; }
    return *this;
  }
  /**
   * \brief Multiply rhs to this series. This is not point wise multiplication,
   * but is convolution.
   *
   * \param rhs
   * \return vec_series&
   */
  vec_series& operator*=(vec_series rhs) {
    // Do FFT only when both of the series is not empty.
    if (!empty() && !rhs.empty()) [[likely]] {
      // Size of LHS divided by 8
      const int lsize_8 = size_8();
      // Size of RHS divided by 8
      const int rsize_8 = rhs.size_8();
      // Size of LHS*RHS divided by 8
      const int nsize_8 = lsize_8 + rsize_8;
      // The number of points to be transformed
      const int points_8 = detail::ceil_pow2(nsize_8);
      // Reise
      c_.resize(points_8);
      rhs.c_.resize(points_8);
      // Do transform
      detail::simd::dif_butterfly(std::span{c_});
      detail::simd::dif_butterfly(std::span{rhs.c_});
      // Element wise multiplication
      for (int i = 0; i < points_8; ++i) { c_[i] *= rhs.c_[i]; }
      // Inverse transform
      detail::simd::dit_butterfly(std::span{c_});
      // inv(points_8 * 8)
      const auto ipoints = mint_type(points_8 * 8).inv();
      const mintx8_type ipoints_x8{_mm256_set1_epi32(ipoints.get())};
      // Normalize
      c_.resize(nsize_8);
      for (auto& f : c_) { f *= ipoints_x8; }
    } else {
      // RHS will be destroyed so we can reuse the resource.
      if (c_.capacity() < rhs.c_.capacity()) { c_.swap(rhs.c_); }
      c_.clear();
    }
    return *this;
  }
  /**
   * \brief Constant multiplication.
   *
   * \param rhs
   * \return vec_series&
   */
  vec_series& operator*=(mint_type rhs) {
    const mintx8_type rhs_x8{_mm256_set1_epi32(rhs.get())};
    for (auto& f : c_) { f *= rhs_x8; }
    return *this;
  }
  mint_type operator[](int i) const {
    if (0 <= i && i < size()) { return c_[i / 8].extract(i % 8); }
    return 0;
  }
  mint_type operator()(mint_type x) const {
    alignas(__m256i) std::array<mint_type, 8> buf;
    mint_type res = 0;
    mint_type fac = 1;
    for (const auto& f : c_) {
      f.store_aligned(buf.data());
      for (int i = 0; i < 8; ++i) {
        res += buf[i] * fac;
        fac *= x;
      }
    }
    return res;
  }
  /**
   * \brief Just return the argument.
   *
   * \param fx
   * \return vec_series
   */
  friend auto operator+(vec_series fx) -> vec_series { return fx; }
  /**
   * \brief Negate given polynomial.
   *
   * \param fx
   * \return vec_series
   */
  friend auto operator-(vec_series fx) -> vec_series {
    for (auto& f : fx.c_) { f = -f; }
    return fx;
  }
  // Binary operaitons
  /**
   * \brief Polynomial addition.
   *
   * \param lhs
   * \param rhs
   * \return vec_series
   */
  friend auto operator+(vec_series lhs, const vec_series& rhs) -> vec_series {
    lhs += rhs;
    return lhs;
  }
  /**
   * \brief Polynomial subtraction.
   *
   * \param lhs
   * \param rhs
   * \return vec_series
   */
  friend auto operator-(vec_series lhs, const vec_series& rhs) -> vec_series {
    lhs -= rhs;
    return lhs;
  }
  /**
   * \brief Polynomial multiplication.
   *
   * \param lhs
   * \param rhs
   * \return vec_series
   */
  friend auto operator*(vec_series lhs, vec_series rhs) -> vec_series {
    lhs *= std::move(rhs);
    return lhs;
  }
  /**
   * \brief
   *
   * \param lhs
   * \param rhs
   * \return vec_series
   */
  friend auto operator*(vec_series lhs, mint_type rhs) -> vec_series {
    lhs *= rhs;
    return lhs;
  }
  /**
   * \brief
   *
   * \param lhs
   * \param rhs
   * \return vec_series
   */
  friend auto operator*(mint_type lhs, vec_series rhs) -> vec_series {
    rhs *= lhs;
    return rhs;
  }
  /**
   * \brief Compute x*F'(x)
   *
   * \return vec_series&
   */
  vec_series& diff_x() {
    mintx8_type fac{_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)};
    const mintx8_type step{_mm256_set1_epi32(mint_type(8).get())};
    fac.ymm = mintx8_type::transform(fac.ymm);
    for (int i = 0; i < size_8(); ++i) {
      c_[i] *= fac;
      fac += step;
    }
    return *this;
  }
  /**
   * \brief Compute integral(F(x))/x
   *
   * \return vec_series&
   */
  vec_series& integrate_x() {
    const auto invs = mint_type::inv_table(size());
    mintx8_type buf;
    for (int i = 0; i < size_8(); ++i) {
      buf.load(&invs[8 * i]);
      c_[i] *= buf;
    }
    return *this;
  }
  /**
   * \brief Compute F(x) * F(x)
   *
   * \return vec_series&
   */
  vec_series& square() {
    if (!empty()) [[likely]] {
      // The size of fx*fx divided by 8
      const int nsize_8 = 2 * size_8();
      // The number of points to be transform
      const int points_8 = detail::ceil_pow2(nsize_8);
      // Resize
      c_.resize(points_8);
      // Do DFT
      detail::simd::dif_butterfly(std::span{c_});
      // Element wise multiplication
      for (auto& f : c_) { f *= f; }
      // Inverse DFT
      detail::simd::dit_butterfly(std::span{c_});
      // Resize & Normalize
      const auto ipoints = mint_type(points_8 * 8).inv();
      const mintx8_type ipoints_x8{_mm256_set1_epi32(ipoints.get())};
      c_.resize(nsize_8);
      for (auto& f : c_) { f *= ipoints_x8; }
    }
    return *this;
  }
  // Common functions
  /**
   * \brief Compute inv(F(x)) in O(NlogN)
   *
   * \param n
   * \return vec_series
   */
  vec_series inv(int n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
    assert((*this)(0) != 0);
#endif
    // The size of result
    const int nsize_8 = std::max(1, detail::ceil_pow2(n) / 8);
    // The result
    std::vector<mintx8_type> q_x(nsize_8);
    // For first 8 terms, we can naively compute it.
    {
      // Reference: https://en.wikipedia.org/wiki/Power_series
      alignas(__m256i) std::array<mint_type, 8> a0;
      alignas(__m256i) std::array<mint_type, 8> b0;
      // Store the first 8 terms into a0.
      c_.front().store_aligned(a0.data());
      // inv(a0[0])
      const auto inv_a0 = a0.front().inv();
      b0.front() = inv_a0;
      for (int i = 1; i < 8; ++i) {
        b0[i] = 0;
        for (int j = 0; j < i; ++j) { b0[i] -= b0[j] * a0[i - j]; }
        b0[i] *= inv_a0;
      }
      q_x.front().load_aligned(b0.data());
    }
    // Step
    const mintx8_type i4_x8{_mm256_set1_epi32(mint_type(4).inv().get())};
    mintx8_type isize_x8{_mm256_set1_epi32(mint_type(64).inv().get())};
    // isize_x8 reflect the inverse of square of the current size.
    for (int k = 1; k < nsize_8; k *= 2) {
      // k refer to the current size divided by 8
      // Update isize_x8
      isize_x8 *= i4_x8;
      // The size we want to compute
      const int cur_k = k;
      const int nxt_k = k * 2;
      // Buffer to store the result of DFT
      std::vector<mintx8_type> buf_fx(nxt_k);
      std::vector<mintx8_type> buf_qx(nxt_k);
      // Q_nxt = Q_prev + Q_prev * (1 - F * Q_prev)
      std::copy_n(c_.begin(), std::min(size_8(), nxt_k), buf_fx.begin());
      std::copy_n(q_x.begin(), k, buf_qx.begin());
      // Do DFT
      detail::simd::dif_butterfly(std::span{buf_fx});
      detail::simd::dif_butterfly(std::span{buf_qx});
      // Element wise multiplication
      for (int i = 0; i < nxt_k; ++i) { buf_fx[i] *= buf_qx[i]; }
      // Inverse DFT
      detail::simd::dit_butterfly(std::span{buf_fx});
      // buf_fx now holds the result of (B*A)-1 mod x^2k
      // We also normalize buf_fx
      std::fill_n(buf_fx.begin(), cur_k, mintx8_type{});
      for (int i = cur_k; i < nxt_k; ++i) { buf_fx[i] *= isize_x8; }
      // Do DFT
      detail::simd::dif_butterfly(std::span{buf_fx});
      // We have already computed DFT for buf_qx so we can reuse the result
      for (int i = 0; i < nxt_k; ++i) { buf_qx[i] *= -buf_fx[i]; }
      // Inverse DFT
      detail::simd::dit_butterfly(std::span{buf_qx});
      // Store the result
      std::copy_n(buf_qx.begin() + cur_k, cur_k, q_x.begin() + cur_k);
    }
    vec_series res{std::move(q_x)};
    res.resize(n);
    return res;
  }
  /**
   * \brief Compute exp(f(x)) mod x^n in O(NlogN).
   *
   * \param n n >= 0 required.
   * \return vec_series
   */
  vec_series exp(int n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
    assert((*this)(0) == 0);
#endif
    /*
      F = *this mod x^k
      Q = exp(F) mod x^k
      R = 1/Q mod x^k
      F_nxt = L + x^k H mod x^2k
      Q_nxt = Q + Q * (F_nxt - log(Q)) mod x^2k
      Q_nxt = Q + Q * I(R * ((x*F_nxt'*Q - x*Q')) mod x^2k
      Q_nxt = Q + Q * I(R * ((x*L*Q-x*Q'))+H*x^k)
    */
    // Note that x*Q' = (x*F')*Q, so we can naively compute the first 8 terms.
    // The size of the result divided by 8
    const int nsize_8 = std::max(1, detail::ceil_pow2(n) / 8);
    // x*F'
    std::vector<mintx8_type> f_x(nsize_8);
    // exp(F(x))
    std::vector<mintx8_type> q_x(nsize_8);
    // inv(Q(x))
    std::vector<mintx8_type> r_x(nsize_8);
    // Inverse table
    const auto invs = mint_type::inv_table(8 * nsize_8);
    // Compute x*F'
    {
      mintx8_type fac{_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)};
      const mintx8_type step{_mm256_set1_epi32(mint_type(8).get())};
      fac.ymm = mintx8_type::transform(fac.ymm);
      for (int i = 0; i < std::min(size_8(), nsize_8); ++i) {
        f_x[i] = c_[i] * fac;
        fac += step;
      }
    }
    // Naively compute first 8 terms of Q and R.
    {
      // q(i) = (1/i)*sum(q(j)*f(i-j) for 0 <= j < i) for i > 0.
      // r(i) = -(1/q(0))*sum(r(j)*q(i-j) for 0<= j < i) for i > 0.
      alignas(__m256i) std::array<mint_type, 8> f0;
      alignas(__m256i) std::array<mint_type, 8> q0;
      alignas(__m256i) std::array<mint_type, 8> r0;
      // Store the first 8 terms into f0
      f_x.front().store_aligned(f0.data());
      // q[0] = 1
      q0.front() = 1;
      // r[0] = 1
      r0.front() = 1;
      for (int i = 1; i < 8; ++i) {
        // First compute q[i]
        q0[i] = 0;
        for (int j = 0; j < i; ++j) { q0[i] += q0[j] * f0[i - j]; }
        q0[i] *= invs[i];
        // Next compute r[i]
        r0[i] = 0;
        for (int j = 0; j < i; ++j) { r0[i] -= r0[j] * q0[i - j]; }
        // r0[i] *= 1
      }
      q_x.front().load_aligned(q0.data());
      r_x.front().load_aligned(r0.data());
    }
    const mintx8_type i2_x8{_mm256_set1_epi32(mint_type(2).inv().get())};
    mintx8_type isize_x8{_mm256_set1_epi32(mint_type(8).inv().get())};
    // Start computing remaining terms
    for (int k = 1; k < nsize_8; k *= 2) {
      // q_x = exp(f_x) mod x^(8k)
      // r_x = inv(q_x) mod x^(8k)
      // Next size
      const int cur_k = k;
      const int nxt_k = k * 2;
      const auto cur_isize_x8 = isize_x8;
      const auto nxt_isize_x8 = isize_x8 * i2_x8;
      // Buffer to store the result of DFT of F
      std::vector<mintx8_type> buf_fx(nxt_k);
      // Buffer to store the result of DFT of Q
      std::vector<mintx8_type> buf_qx(nxt_k);
      // Buffer to store the result of DFT of R
      std::vector<mintx8_type> buf_rx(nxt_k);
      // Copy data
      std::copy_n(f_x.begin(), cur_k, buf_fx.begin());
      std::copy_n(q_x.begin(), cur_k, buf_qx.begin());
      std::copy_n(r_x.begin(), cur_k, buf_rx.begin());
      // Transform
      detail::simd::dif_butterfly(
          std::span{buf_fx.begin(), buf_fx.begin() + cur_k});
      detail::simd::dif_butterfly(
          std::span{buf_qx.begin(), buf_qx.begin() + nxt_k});
      detail::simd::dif_butterfly(
          std::span{buf_rx.begin(), buf_rx.begin() + nxt_k});
      // buf_fx = (x * L' * Q) mod x^k - 1.
      // This follows that ([x*L'*Q mod x^k - 1] - x*Q')*x^k = (x*L'*Q - x*Q')
      // because x*L'*Q = x*Q' mod x^k
      for (int i = 0; i < cur_k; ++i) { buf_fx[i] *= buf_qx[i]; }
      // Inverse transform
      detail::simd::dit_butterfly(
          std::span{buf_fx.begin(), buf_fx.begin() + cur_k});
      // Now buf_fx = (x*L'*Q - x*Q')/x^k
      {
        mintx8_type fac{_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)};
        const mintx8_type step{_mm256_set1_epi32(mint_type(8).get())};
        fac.ymm = mintx8_type::transform(fac.ymm);
        for (int i = 0; i < cur_k; ++i) {
          buf_fx[i] *= cur_isize_x8;
          buf_fx[i] -= fac * q_x[i];
          fac += step;
        }
      }
      // Now buf_fx = (x*L'*Q - x*Q')
      std::copy_n(buf_fx.begin(), cur_k, buf_fx.begin() + cur_k);
      std::fill_n(buf_fx.begin(), cur_k, mintx8_type{});
      // Multiply by R
      detail::simd::dif_butterfly(std::span{buf_fx});
      // Element wise mul
      for (int i = 0; i < nxt_k; ++i) { buf_fx[i] *= buf_rx[i]; }
      // Inverse transform
      detail::simd::dit_butterfly(std::span{buf_fx});
      // Now buf_fx = R*(x*L'*Q-x*Q')
      std::fill_n(buf_fx.begin(), cur_k, mintx8_type{});
      // Now buf_fx = I(R*(x*L'*Q - x*Q')+k*H'*x^k)
      {
        mintx8_type buf;
        for (int i = cur_k; i < nxt_k; ++i) {
          // Normalize
          buf_fx[i] *= nxt_isize_x8;
          // Add k*H'*x^k
          buf_fx[i] += f_x[i];
          buf.load(&invs[8 * i]);
          // Integrate
          buf_fx[i] *= buf;
        }
      }
      // Multiply by Q
      detail::simd::dif_butterfly(std::span{buf_fx});
      // The result is stored in buf_qx
      for (int i = 0; i < nxt_k; ++i) { buf_qx[i] *= buf_fx[i]; }
      // Inverse transform
      detail::simd::dit_butterfly(std::span{buf_qx});
      // Normalize the result
      for (int i = cur_k; i < nxt_k; ++i) { buf_qx[i] *= nxt_isize_x8; }
      // Store the result
      std::copy_n(buf_qx.begin() + cur_k, cur_k, q_x.begin() + cur_k);
      if (nxt_k != nsize_8) {
        // Update r_x
        std::copy_n(q_x.begin(), cur_k, buf_qx.begin());
        detail::simd::dif_butterfly(std::span{buf_qx});
        // Element wise multiply
        for (int i = 0; i < nxt_k; ++i) { buf_qx[i] *= buf_rx[i]; }
        // Inverse transform
        detail::simd::dit_butterfly(std::span{buf_qx});
        // Normalize
        std::fill_n(buf_qx.begin(), cur_k, mintx8_type{});
        for (int i = cur_k; i < nxt_k; ++i) { buf_qx[i] *= nxt_isize_x8; }
        // Transform
        detail::simd::dif_butterfly(std::span{buf_qx});
        // Element wise mul
        // The result is stored in buf_rx
        for (int i = 0; i < nxt_k; ++i) { buf_rx[i] *= -buf_qx[i]; }
        // Inverse transform
        detail::simd::dit_butterfly(std::span{buf_rx});
        // Normalize
        for (int i = cur_k; i < nxt_k; ++i) { buf_rx[i] *= nxt_isize_x8; }
        // Store the result
        std::copy_n(buf_rx.begin() + cur_k, cur_k, r_x.begin() + cur_k);
        // Update isize
        isize_x8 = nxt_isize_x8;
      }
    }
    vec_series res{std::move(q_x)};
    res.resize(n);
    return res;
  }
  /**
   * \brief Compute log(F(x)) mod x^n in O(NlogN).
   *
   * \param n n >= 0 required.
   * \return vec_series
   */
  vec_series log(int n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
    assert((*this)(0) != 0);
#endif
    auto res{*this};
    res.resize(n);
    // log(f(x)) = integral(diff(f(x))/f(x))
    // log(f(x)) = integral(x*diff(f(x))/f(x))/x
    {
      mintx8_type fac{_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)};
      const mintx8_type step{_mm256_set1_epi32(mint_type(8).get())};
      fac.ymm = mintx8_type::transform(fac.ymm);
      for (int i = 0; i < res.size_8(); ++i) {
        res.c_[i] *= fac;
        fac += step;
      }
    }
    res *= inv(n);
    res.resize(n);
    {
      const auto invs = mint_type::inv_table(res.size());
      mintx8_type buf;
      for (int i = 0; i < res.size_8(); ++i) {
        buf.load(&invs[8 * i]);
        res.c_[i] *= buf;
      }
    }
    return res;
  }
  /**
   * \brief Resize the series.
   *
   * \param n
   */
  void resize(int n) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    detail::assume(n >= 0);
    const int nsize_8 = (n + 7) / 8;
    c_.resize(nsize_8);
    if (n % 8 != 0) { c_.back().cutoff(n % 8); }
  }
  /**
   * \brief Return the size of this series.
   *
   * \return int
   */
  int size() const { return size_8() * 8; }
  /**
   * \brief Check if this series is empty or not.
   *
   * \return true
   * \return false
   */
  bool empty() const { return c_.empty(); }
  /**
   * \brief Retrieve the coefficients data. Return the vector with size n where
   * each element at index i is equal to coef(i).
   *
   * \return std::vector<mint_type>
   */
  std::vector<mint_type> to_vector() const {
    std::vector<mint_type> res(size());
    // Temporary buffer to store each packed data
    alignas(__m256i) std::array<mint_type, 8> buf;
    for (int i = 0; i < size_8(); ++i) {
      c_[i].store_aligned(buf.data());
      for (int j = 0; j < 8; ++j) { res[i * 8 + j] = buf[j]; }
    }
    return res;
  }
  /**
   * \brief Return the internal data.
   *
   * \return std::vector<mintx8_type>
   */
  std::vector<mintx8_type> data() const& { return c_; }
  /**
   * \brief Return the internal data.
   *
   * \return std::vector<mintx8_type>
   */
  std::vector<mintx8_type> data() && { return std::move(c_); }
private:
  int size_8() const { return static_cast<int>(c_.size()); }
  // Coefficients data.
  std::vector<mintx8_type> c_;
};
}  // namespace algo