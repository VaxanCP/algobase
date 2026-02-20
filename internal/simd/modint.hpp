#pragma once

// Low level simd modint.

#include <immintrin.h>

#include <array>
#include <bit>
#include <vector>

#include "../../class/montgomery.hpp"
#include "../modular.hpp"
#include "./helpers.hpp"

/*
@internal/simd/helpers.hpp
@internal/modular.hpp
@class/montgomery.hpp
*/

// makecode
namespace algo::detail::simd {
/**
 * \brief Low level vectorized montgomery modint implementation. You should use
 * the class only when you know what you are doing.
 *
 * \tparam MOD
 */
template <uint32_t MOD = 998244353>
struct montgomery_modintx8 {
  using mint_type = montgomery_modint<MOD>;
  static_assert(MOD % 2 == 1, "Modulus must be odd");
  __m256i ymm;
  // No constructors
  montgomery_modintx8& operator+=(montgomery_modintx8 rhs) {
    ymm = add_mod(ymm, rhs.ymm);
    return *this;
  }
  montgomery_modintx8& operator-=(montgomery_modintx8 rhs) {
    ymm = sub_mod(ymm, rhs.ymm);
    return *this;
  }
  montgomery_modintx8& operator*=(montgomery_modintx8 rhs) {
    ymm = mul_mod(ymm, rhs.ymm);
    return *this;
  }
  friend auto operator+(montgomery_modintx8 d) -> montgomery_modintx8 {
    return d;
  }
  friend auto operator-(montgomery_modintx8 d) -> montgomery_modintx8 {
    d.ymm = mask_negate<0b11111111>(d.ymm);
    d.ymm = normalize(d.ymm);
    return d;
  }
  friend auto operator+(montgomery_modintx8 lhs, montgomery_modintx8 rhs)
      -> montgomery_modintx8 {
    lhs += rhs;
    return lhs;
  }
  friend auto operator-(montgomery_modintx8 lhs, montgomery_modintx8 rhs)
      -> montgomery_modintx8 {
    lhs -= rhs;
    return lhs;
  }
  friend auto operator*(montgomery_modintx8 lhs, montgomery_modintx8 rhs)
      -> montgomery_modintx8 {
    lhs *= rhs;
    return lhs;
  }
  /**
   * \brief Extract i-th value of the vector.
   *
   * \param i 0 <= i < 8 must hold
   * \return mint_type
   */
  mint_type extract(int i) const {
    alignas(__m256i) mint_type buf[8];
    store_aligned(buf);
    return buf[i];
  }
  /**
   * \brief Change i-th element of the vector.
   *
   * \param i 0 <= i < 8 must hold.
   * \param x
   */
  void insert(int i, mint_type x) {
    alignas(__m256i) mint_type buf[8];
    store_aligned(buf);
    buf[i] = x;
    load_aligned(buf);
  }
  /**
   * \brief Cut off vector to n elements. The last 8-n elements are set to zero.
   *
   * \param n 0 <= n <= 8 must hold.
   */
  void cutoff(int n) {
    // Set the last 8-n elements to zero
    const int32_t maskl[16]{-1, -1, -1, -1, -1, -1, -1, -1,
                            0,  0,  0,  0,  0,  0,  0,  0};
    const auto xmm =
        _mm256_loadu_si256(reinterpret_cast<const __m256i_u*>(maskl + 8 - n));
    ymm = _mm256_and_si256(ymm, xmm);
  }
  /**
   * \brief Load and create new montgomery_modintx8. The address need not to be
   * aligned.
   *
   * \param mem
   * \return montgomery_modintx8
   */
  void load(const void* mem) {
    ymm = _mm256_loadu_si256(reinterpret_cast<const __m256i_u*>(mem));
  }
  /**
   * \brief Load data from given address. The address must be aligned.
   *
   * \param mem
   * \return montgomery_modintx8
   */
  void load_aligned(const void* mem) {
    ymm = _mm256_load_si256(reinterpret_cast<const __m256i*>(mem));
  }
  /**
   * \brief Store this->ymm into the provided address. The address need not to
   * be aligned.
   *
   * \param mem
   */
  void store(void* mem) const {
    _mm256_storeu_si256(reinterpret_cast<__m256i_u*>(mem), ymm);
  }
  /**
   * \brief Store this->ymm into the provided address. The address have to be
   * aligned.
   *
   * \param mem
   */
  void store_aligned(void* mem) const {
    _mm256_store_si256(reinterpret_cast<__m256i*>(mem), ymm);
  }
  /**
   * \brief A[i] += (A[i] < 0) ? MOD : 0;
   *
   * \param A
   * \return __m256i
   */
  static __m256i normalize(__m256i A) {
    const auto md = _mm256_set1_epi32(MOD);
    const auto mask = _mm256_srai_epi32(A, 31);
    return _mm256_add_epi32(A, _mm256_and_si256(md, mask));
  }
  /**
   * \brief Add two vectors modulo MOD
   *
   * \param A
   * \param B
   * \return __m256i
   */
  static __m256i add_mod(__m256i A, __m256i B) {
    const auto md = _mm256_set1_epi32(MOD);
    // A + B
    const auto d = _mm256_add_epi32(A, B);
    // A + B - M
    const auto r = _mm256_sub_epi32(d, md);
    return normalize(r);
  }
  /**
   * \brief Subtract one vector from the other modulo MOD
   *
   * \param A
   * \param B
   * \return __m256i
   */
  static __m256i sub_mod(__m256i A, __m256i B) {
    // A - B
    const auto d = _mm256_sub_epi32(A, B);
    return normalize(d);
  }
  static __m256i mul_mod(__m256i A, __m256i B) {
    // MOD^(-1) mod 2^32
    // Montgomery multiplication. For detailed information, see here:
    // https://cp-algorithms.com/algebra/montgomery_multiplication.html
    constexpr uint32_t Nd = mint_type::Nd;
    const auto md = _mm256_set1_epi32(MOD);
    const auto nd = _mm256_set1_epi32(Nd);
    const auto prod0246 = _mm256_mul_epu32(A, B);
    A = _mm256_srli_si256(A, 4);
    B = _mm256_srli_si256(B, 4);
    const auto prod1357 = _mm256_mul_epu32(A, B);
    const auto qlo = _mm256_mul_epu32(prod0246, nd);
    const auto qhi = _mm256_mul_epu32(prod1357, nd);
    const auto tmplo = _mm256_mul_epu32(qlo, md);
    const auto tmphi = _mm256_mul_epu32(qhi, md);
    const auto alo = _mm256_sub_epi32(prod0246, tmplo);
    const auto ahi = _mm256_sub_epi32(prod1357, tmphi);
    const auto tmpd = shuffle_x32<_MM_SHUFFLE(3, 1, 3, 1)>(alo, ahi);
    const auto d = permute_x32<_MM_SHUFFLE(3, 1, 2, 0)>(tmpd);
    return normalize(d);
  }
  /**
   * \brief Multiply A with (2^32)^(-1) mod MOD
   *
   * \param A
   * \return __m256i
   */
  static __m256i reduce(__m256i A) {
    constexpr uint32_t Nd = mint_type::Nd;
    const auto md = _mm256_set1_epi32(MOD);
    const auto nd = _mm256_set1_epi32(Nd);
    const auto qlo = _mm256_mul_epu32(A, nd);
    A = _mm256_srli_si256(A, 4);
    const auto qhi = _mm256_mul_epu32(A, nd);
    const auto prodlo = _mm256_mul_epu32(qlo, md);
    const auto prodhi = _mm256_mul_epu32(qhi, md);
    const auto tmplo = _mm256_sub_epi32(_mm256_setzero_si256(), prodlo);
    const auto tmphi = _mm256_sub_epi32(_mm256_setzero_si256(), prodhi);
    const auto tmpd = shuffle_x32<_MM_SHUFFLE(3, 1, 3, 1)>(tmplo, tmphi);
    const auto d = permute_x32<_MM_SHUFFLE(3, 1, 2, 0)>(tmpd);
    return normalize(d);
  }
  /**
   * \brief Multiply A with 2^64 in montgomery space.
   *
   * \param A
   * \return __m256i
   */
  static __m256i transform(__m256i A) {
    constexpr uint32_t R2 = mint_type::R2;
    const auto rsq = _mm256_set1_epi32(R2);
    return mul_mod(A, rsq);
  }
};
}  // namespace algo::detail::simd