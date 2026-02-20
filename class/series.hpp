#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <span>
#include <vector>

#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
#include "../internal/dif.hpp"
/*
@internal/base/bit-base.hpp
@internal/base/typing.hpp
@internal/base/def.hpp
@internal/dif.hpp
*/

// makecode
namespace algo {
/**
 * \brief Fast formal power series implementation. This class is specialized
 * implementation for 'NTT friendly modulus' like a 998244353. In formal power
 * series, we don't care about the degree of polynomial at all. If you want to
 * reduce the memory usage, call c.resize(). We assume that during some
 * operations, the size of 'c' doesn't exceed ModNum::mod() and specifically
 * for multiplication, the size of 'c' doesn't exceed
 * 2^(count_tz(ModNum::mod()-1)).
 *
 * \tparam ModNum
 */
template <detail::modular_integer ModNum>
  requires(ModNum::is_prime_mod())
struct series {
  // Coefficients data.
  std::vector<ModNum> c;
  /**
   * \brief Get the i th coefficient of the data if exist, or return 0
   * otherwise.
   *
   * \param i
   * \return ModNum
   */
  ModNum operator[](int i) const {
    if (0 <= i && i < size()) { return c[i]; }
    return 0;
  }
  /**
   * \brief Evaluate f(x)
   *
   * \param x
   * \return ModNum
   */
  ModNum operator()(ModNum x) const {
    ModNum res = 0;
    ModNum fac = 1;
    for (const auto v : c) {
      res += fac * v;
      fac *= x;
    }
    return res;
  }
  /**
   * \brief Add two polynomial
   *
   * \param rhs
   * \return fast_fps&
   */
  series& operator+=(const series& rhs) {
    // Resize to the larger one
    if (c.size() < rhs.c.size()) { c.resize(rhs.c.size()); }
    // Do addition
    for (int i = 0; i < rhs.size(); ++i) { c[i] += rhs.c[i]; }
    return *this;
  }
  /**
   * \brief Subtract lhs from rhs.
   *
   * \param rhs
   * \return fast_fps&
   */
  series& operator-=(const series& rhs) {
    // Resize to the larger one
    if (c.size() < rhs.c.size()) { c.resize(rhs.c.size()); }
    // Do subtraction
    for (int i = 0; i < rhs.size(); ++i) { c[i] -= rhs.c[i]; }
    return *this;
  }
  /**
   * \brief Multiply two polynomials.
   *
   * \param rhs
   * \return fast_fps&
   */
  series& operator*=(series rhs) {
    // Check if both of them are not empty
    if (!empty() && !rhs.empty()) [[likely]] {
      // The size of lhs
      const int size_l = size();
      // The size of rhs
      const int size_r = rhs.size();
      // Size of (lhs*rhs)
      const int size_lr = size_l + size_r - 1;
      // The number of points to be transformed
      const int num_points = detail::ceil_pow2(size_lr);
      // Resize
      c.resize(num_points);
      rhs.c.resize(num_points);
      // Do DFT
      detail::dif_butterfly(std::span{c});
      detail::dif_butterfly(std::span{rhs.c});
      // Element wise multiplication
      for (int i = 0; i < num_points; ++i) { c[i] *= rhs.c[i]; }
      // Do inverse transform
      detail::dit_butterfly(std::span{c});
      // Resize
      c.resize(size_lr);
      // Invers of num_points
      const auto inum_points = ModNum(num_points).inv();
      // Normalize the result
      for (auto& f : c) { f *= inum_points; }
    } else {
      // Rhs will be destroyed after the call, so we can use swap the resource
      // if needed.
      if (c.capacity() < rhs.c.capacity()) { c.swap(rhs.c); }
      // Since the result will be zero
      c.clear();
    }
    return *this;
  }
  /**
   * \brief Multiply *this by a constant.
   *
   * \param rhs
   * \return fast_fps&
   */
  series& operator*=(ModNum rhs) {
    // Constant multiplication
    for (auto& f : c) { f *= rhs; }
    return *this;
  }
  /**
   * \brief Multiply this series by x^n
   *
   * \param n
   * \return series&
   */
  series& operator<<=(int n) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    c.resize(size() + n);
    std::shift_right(c.begin(), c.end(), n);
    std::fill(c.begin(), c.begin() + n, 0);
    return *this;
  }
  /**
   * \brief
   *
   * \param n
   * \return series&
   */
  series& operator>>=(int n) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    if (size() > n) [[likely]] {
      std::shift_left(c.begin(), c.end(), n);
      c.resize(size() - n);
    } else {
      c.clear();
    }
    return *this;
  }
  /**
   * \brief Unary plus operator. Just return the argument.
   *
   * \param fx
   * \return fast_fps
   */
  friend series operator+(series fx) { return fx; }
  /**
   * \brief Unary minus operator. Negate all elements of the given series.
   *
   * \param fx
   * \return fast_fps
   */
  friend series operator-(series fx) {
    for (auto& f : fx.c) { f = -f; }
    return fx;
  }
  /**
   * \brief Add two polynomials.
   *
   * \param lhs
   * \param rhs
   * \return fast_fps
   */
  friend series operator+(series lhs, const series& rhs) {
    lhs += rhs;
    return lhs;
  }
  /**
   * \brief Perform subtraction between two polynomials.
   *
   * \param lhs
   * \param rhs
   * \return fast_fps
   */
  friend series operator-(series lhs, const series& rhs) {
    lhs -= rhs;
    return lhs;
  }
  /**
   * \brief Multiply two polynomials.
   *
   * \param lhs
   * \param rhs
   * \return fast_fps
   */
  friend series operator*(series lhs, series rhs) {
    lhs *= std::move(rhs);
    return lhs;
  }
  /**
   * \brief Multiply LHS by a constant RHS.
   *
   * \param lhs
   * \param rhs
   * \return fast_fps
   */
  friend series operator*(series lhs, ModNum rhs) {
    lhs *= rhs;
    return lhs;
  }
  /**
   * \brief Multiply RHS by a constant LHS.
   *
   * \param lhs
   * \param rhs
   * \return series
   */
  friend series operator*(ModNum lhs, series rhs) {
    rhs *= lhs;
    return rhs;
  }
  /**
   * \brief Count the number of leading zeroes.
   *
   * \return int
   */
  int leading_zeros() const {
    const auto pos =
        std::find_if(c.rbegin(), c.rend(), [](ModNum x) { return x != 0; });
    return static_cast<int>(c.end() - pos.base());
  }
  /**
   * \brief Count the number of trailing zeroes.
   *
   * \return int
   */
  int trailing_zeros() const {
#if !defined(NDEBUG)
    assert((*this)(0) != 0);
#endif
    const auto pos =
        std::find_if(c.begin(), c.end(), [](ModNum x) { return x != 0; });
    return static_cast<int>(pos - c.begin());
  }
  /**
   * \brief Resize
   *
   * \param n
   */
  void resize(int n) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    c.resize(n);
  }
  /**
   * \brief The current size of FPS
   *
   * \return int
   */
  int size() const { return static_cast<int>(c.size()); }
  /**
   * \brief Check if the series is empty (i.e., 0) or not.
   *
   * \return true
   * \return false
   */
  bool empty() const { return c.empty(); }
};
/**
 * \brief
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> truncate(series<ModNum> fx) {
  fx.resize(fx.size() - fx.leading_zeros());
  return fx;
}
/**
 * \brief Compute f(1/x)*x^deg(f).
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> reverse(series<ModNum> fx) {
  fx = truncate(std::move(fx));
  std::reverse(fx.c.begin(), fx.c.end());
  return fx;
}
/**
 * \brief Compute x*diff(fx).
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> diff_x(series<ModNum> fx) {
  for (int i = 0; i < fx.size(); ++i) { fx.c[i] *= ModNum(i); }
  return fx;
}
/**
 * \brief Differentiate given series.
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> diff(series<ModNum> fx) {
  if (!fx.empty()) [[likely]] {
    // c[i] *= i
    fx = diff_x(std::move(fx));
    // Shift by 1
    std::shift_left(fx.c.begin(), fx.c.end(), 1);
    // Then pop back
    fx.c.pop_back();
  }
  return fx;
}
/**
 * \brief Compute integral(fx) / x.
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> integrate_x(series<ModNum> fx) {
  // c[i] *= 1/i
  const auto invs = ModNum::inv_table(fx.size());
  for (int i = 0; i < fx.size(); ++i) { fx.c[i] *= invs[i]; }
  return fx;
}
/**
 * \brief Integrate given series.
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> integrate(series<ModNum> fx) {
  // integrate_x(fx) * x
  fx.c.emplace_back(0);
  // Shift right by 1
  std::shift_right(fx.c.begin(), fx.c.end(), 1);
  // fx.c[0] = 0
  fx.c.front() = 0;
  return integrate_x(std::move(fx));
}
/**
 * \brief Compute f*f
 *
 * \tparam ModNum
 * \param fx
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> square(series<ModNum> fx) {
  if (!fx.empty()) [[likely]] {
    const int nsize = 2 * fx.size() - 1;
    const int points = detail::ceil_pow2(nsize);
    // Resize
    fx.resize(points);
    // Transform
    detail::dif_butterfly(std::span{fx.c});
    // Mul
    for (auto& f : fx.c) { f *= f; }
    // Inverse transform
    detail::dit_butterfly(std::span{fx.c});
    const auto ipoints = ModNum(points).inv();
    // Normalize
    fx.resize(nsize);
    for (auto& f : fx.c) { f *= ipoints; }
  }
  return fx;
}
/**
 * \brief Compute inverse series of given series up to n coefficients.
 *
 * \tparam ModNum
 * \param fx
 * \param n
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> inv(const series<ModNum>& fx, int n) {
#if !defined(NDEBUG)
  assert(n >= 0);
  assert(fx(0) != 0);
#endif
  // For algorithm itself, please refer to
  // here:https://cp-algorithms.com/algebra/polynomial.html
  /*
    Note:
    When we convolute two polynomials A and B with degree less than n using
    dft of order n, we get the A * B mod x^n-1. If deg(A) < n/2 and deg(B) <
    n, then the result of cyclic convolution between A and B is such a
    polynomial that the coefficients n/2,n/2+1,...n-1 is indeed the
    coefficients of multiplication of A and B.
    // For father reference, see here:
    https://math.mit.edu/classes/18.783/2021/LectureNotes3.pdf
    // A notable consequence of this fact is that for two polynomials A and B,
    DFT(A * B mod x^n-1) is equal to DFT(A,n)*DFT(B,n).
  */
  // Resulting size
  const int nsize = detail::ceil_pow2(n);
  // Suppose we have polynomial A with degree less than k and polynomial B
  // with degree less than l, and we convolute them with size n, then at most
  // first k+l-n-1 terms of result are invalid.
  // The result
  std::vector<ModNum> q_x(nsize);
  // 2^(-1)
  const auto i4 = ModNum(4).inv();
  // The inverse of current size
  ModNum isize = 1;
  // The result of 2^0
  q_x[0] = fx.c.front().inv();
  for (int k = 1; k < nsize; k *= 2) {
    // k is the number of coefficients found so far.
    // Update isize
    isize *= i4;
    // The size we want to compute
    const int cur_k = k;
    const int nxt_k = k * 2;
    // Temporary buffer to store the result of DFT
    // Buffer to store the result of DFT of *this
    std::vector<ModNum> buf_fx(nxt_k);
    std::vector<ModNum> buf_qx(nxt_k);
    // Copy data
    std::copy_n(fx.c.begin(), std::min(fx.size(), nxt_k), buf_fx.begin());
    std::copy_n(q_x.begin(), k, buf_qx.begin());
    // Do DFT
    detail::dif_butterfly(std::span{buf_fx});
    detail::dif_butterfly(std::span{buf_qx});
    // Element wise multiplication
    for (int i = 0; i < nxt_k; ++i) { buf_fx[i] *= buf_qx[i]; }
    // Inverse DFT
    detail::dit_butterfly(std::span{buf_fx});
    // buf_fx now holds the result of (B*A)-1 mod x^2k
    // We also normalize buf_fx
    std::fill_n(buf_fx.begin(), cur_k, 0);
    for (int i = cur_k; i < nxt_k; ++i) { buf_fx[i] *= isize; }
    // Do DFT
    detail::dif_butterfly(std::span{buf_fx});
    // We have already computed DFT for buf_qx, so we can reuse the result
    for (int i = 0; i < nxt_k; ++i) { buf_fx[i] *= -buf_qx[i]; }
    // Inverse DFT
    detail::dit_butterfly(std::span{buf_fx});
    // Store the result
    std::copy_n(buf_fx.begin() + cur_k, cur_k, q_x.begin() + cur_k);
  }
  q_x.resize(n);
  return {std::move(q_x)};
}
/**
 * \brief Compute exp(F) mod x^n.
 *
 * \tparam ModNum
 * \param fx
 * \param n
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> exp(series<ModNum> fx, int n) {
#if !defined(NDEBUG)
  assert(n >= 0);
  assert(fx(0) == 0);
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
  const int nsize = detail::ceil_pow2(n);
  // exp(F(x))
  std::vector<ModNum> q_x(nsize / 1);
  // inv(Q(x))
  std::vector<ModNum> r_x(nsize / 2);
  // exp(F(x)) mod x^1
  q_x.front() = 1;
  // inv(q_x) mod x^1
  if (!r_x.empty()) { r_x.front() = 1; }
  // x*F'
  auto f_x = diff_x(std::move(fx)).c;
  f_x.resize(nsize);
  // Inverse table
  const auto invs = ModNum::inv_table(nsize);
  const auto i2 = ModNum(2).inv();
  ModNum isize = 1;
  for (int k = 1; k < nsize; k *= 2) {
    // q_x = exp(f_x) mod x^k
    // r_x = inv(q_x) mod x^k
    // Next size
    const int cur_k = k;
    const int nxt_k = k * 2;
    // Buffer to store the result of DFT of F
    std::vector<ModNum> buf_fx(nxt_k);
    // Buffer to store the result of DFT of Q
    std::vector<ModNum> buf_qx(nxt_k);
    // Buffer to store the result of DFT of R
    std::vector<ModNum> buf_rx(nxt_k);
    // Copy data
    std::copy_n(f_x.begin(), cur_k, buf_fx.begin());
    std::copy_n(q_x.begin(), cur_k, buf_qx.begin());
    std::copy_n(r_x.begin(), cur_k, buf_rx.begin());
    // Transform
    detail::dif_butterfly(std::span{buf_fx.begin(), buf_fx.begin() + cur_k});
    detail::dif_butterfly(std::span{buf_qx.begin(), buf_qx.begin() + nxt_k});
    detail::dif_butterfly(std::span{buf_rx.begin(), buf_rx.begin() + nxt_k});
    // buf_fx = (x * L' * Q) mod x^k - 1.
    // This follows that ([x*L'*Q mod x^k - 1] - x*Q')*x^k = (x*L'*Q - x*Q')
    // because x*L'*Q = x*Q' mod x^k
    for (int i = 0; i < cur_k; ++i) { buf_fx[i] *= buf_qx[i]; }
    // Inverse transform
    detail::dit_butterfly(std::span{buf_fx.begin(), buf_fx.begin() + cur_k});
    // Now buf_fx = (x*L'*Q - x*Q')/x^k
    for (int i = 0; i < cur_k; ++i) { buf_fx[i] *= isize; }
    for (int i = 0; i < cur_k; ++i) { buf_fx[i] -= ModNum(i) * q_x[i]; }
    // Now buf_fx = (x*L'*Q - x*Q')
    std::copy_n(buf_fx.begin(), cur_k, buf_fx.begin() + cur_k);
    std::fill_n(buf_fx.begin(), cur_k, 0);
    // Multiply by R
    detail::dif_butterfly(std::span{buf_fx});
    // Element wise mul
    for (int i = 0; i < nxt_k; ++i) { buf_fx[i] *= buf_rx[i]; }
    // Inverse transform
    detail::dit_butterfly(std::span{buf_fx});
    // Now buf_fx = R*(x*L'*Q - x*Q')
    std::fill_n(buf_fx.begin(), cur_k, 0);
    // Update isize
    isize *= i2;
    // Normalize
    for (int i = cur_k; i < nxt_k; ++i) { buf_fx[i] *= isize; }
    // Add x*H'*x^k
    // Now buf_fx = R*(x*L'*Q - x*Q')+k*H*x^k
    for (int i = cur_k; i < nxt_k; ++i) { buf_fx[i] += f_x[i]; }
    // Integrate
    // Now buf_fx = I(R*(x*L'*Q - x*Q') + k*H*x^k)
    for (int i = cur_k; i < nxt_k; ++i) { buf_fx[i] *= invs[i]; }
    // Multiply by Q
    detail::dif_butterfly(std::span{buf_fx});
    // The result is stored in buf_qx
    for (int i = 0; i < nxt_k; ++i) { buf_qx[i] *= buf_fx[i]; }
    // Inverse transform.
    detail::dit_butterfly(std::span{buf_qx});
    // Normalize the result
    for (int i = cur_k; i < nxt_k; ++i) { buf_qx[i] *= isize; }
    // Store the result
    std::copy_n(buf_qx.begin() + cur_k, cur_k, q_x.begin() + cur_k);
    if (nxt_k != nsize) {
      // Update r_x
      std::copy_n(q_x.begin(), cur_k, buf_qx.begin());
      detail::dif_butterfly(std::span{buf_qx});
      // Element wise multiply
      for (int i = 0; i < nxt_k; ++i) { buf_qx[i] *= buf_rx[i]; }
      // Inverse transform
      detail::dit_butterfly(std::span{buf_qx});
      // Normalize
      std::fill_n(buf_qx.begin(), cur_k, 0);
      for (int i = cur_k; i < nxt_k; ++i) { buf_qx[i] *= isize; }
      // Transform
      detail::dif_butterfly(std::span{buf_qx});
      // Element wise mul
      // The result is stored in buf_rx
      for (int i = 0; i < nxt_k; ++i) { buf_rx[i] *= -buf_qx[i]; }
      // Inverse transform
      detail::dit_butterfly(std::span{buf_rx});
      // Normalize
      for (int i = cur_k; i < nxt_k; ++i) { buf_rx[i] *= isize; }
      // Store the result
      std::copy_n(buf_rx.begin() + cur_k, cur_k, r_x.begin() + cur_k);
      // Update isize
    }
  }
  q_x.resize(n);
  return {std::move(q_x)};
}
/**
 * \brief Compute log(f(x)) mod x^n.
 *
 * \tparam ModNum
 * \param fx
 * \param n
 * \return series<ModNum>
 */
template <typename ModNum>
series<ModNum> log(series<ModNum> fx, int n) {
#if !defined(NDEBUG)
  assert(n >= 0);
  assert(fx(0) != 0);
#endif
  fx.resize(n);
  auto qx = inv(fx, n);
  qx *= diff_x(std::move(fx));
  qx.resize(n);
  return integrate_x(std::move(qx));
}
}  // namespace algo