#pragma once

#include <immintrin.h>

#include "../class/montgomery.hpp"

/*
@class/montgomery.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Mix and shuffle two vectors. The operation is equivalent to:
 * res[0] = A[0 + N[0:1]]
 * res[1] = A[0 + N[2:3]]
 * res[2] = B[0 + N[4:5]]
 * res[3] = B[0 + N[6:7]]
 * res[4] = A[4 + N[0:1]]
 * res[5] = A[4 + N[2:3]]
 * res[6] = B[4 + N[4:5]]
 * res[7] = B[4 + N[6:7]]
 *
 * \tparam N
 * \param A
 * \param B
 * \return __m256i
 */
template <int N>
static __m256i shuffle256x32(__m256i A, __m256i B) {
  return _mm256_castps_si256(
      _mm256_shuffle_ps(_mm256_castsi256_ps(A), _mm256_castsi256_ps(B), N));
}
/**
 * \brief mix and shuffle given two vectors. The operation is equivalent to:
 * res[0] = A[N[0:1]]
 * res[1] = A[N[2:3]]
 * res[2] = B[N[4:5]]
 * res[3] = B[N[6:7]]
 *
 * \tparam N
 * \param A
 * \param B
 * \return __m128i
 */
template <int N>
static __m128i shuffle128x32(__m128i A, __m128i B) {
  return _mm_castps_si128(
      _mm_shuffle_ps(_mm_castsi128_ps(A), _mm_castsi128_ps(B), N));
}
// DEFINE SELECT4(A, control)
// 	CASE(control[1:0]) OF
// 	0:	tmp[63:0] := A[63:0]
// 	1:	tmp[63:0] := A[127:64]
// 	2:	tmp[63:0] := A[191:128]
// 	3:	tmp[63:0] := A[255:192]
// 	ESAC
// 	RETURN tmp[63:0]
// }
// dst[63:0] := SELECT4(A[255:0], imm8[1:0])
// dst[127:64] := SELECT4(A[255:0], imm8[3:2])
// dst[191:128] := SELECT4(A[255:0], imm8[5:4])
// dst[255:192] := SELECT4(A[255:0], imm8[7:6])
// dst[MAX:256] := 0
template <int N>
static __m256i shuffle256x64(__m256i A) {
  return _mm256_castpd_si256(_mm256_permute4x64_pd(_mm256_castsi256_pd(A), N));
}
// DEFINE SELECT4(src1, src2, control) {
// 	CASE(control[1:0]) OF
// 	0:	tmp[127:0] := src1[127:0]
// 	1:	tmp[127:0] := src1[255:128]
// 	2:	tmp[127:0] := src2[127:0]
// 	3:	tmp[127:0] := src2[255:128]
// 	ESAC
// 	IF control[3]
// 		tmp[127:0] := 0
// 	FI
// 	RETURN tmp[127:0]
// }
// dst[127:0] := SELECT4(a[255:0], b[255:0], imm8[3:0])
// dst[255:128] := SELECT4(a[255:0], b[255:0], imm8[7:4])
// dst[MAX:256] := 0
template <int N>
static __m256i shuffle256x128(__m256i A, __m256i B) {
  return _mm256_castps_si256(_mm256_permute2f128_ps(_mm256_castsi256_ps(A),
                                                    _mm256_castsi256_ps(B), N));
}
/**
 * \brief Permute given vector. The operation is equivalent to:
 * res[0] = A[N[0:1]]
 * res[1] = A[N[2:3]]
 * res[2] = A[N[4:5]]
 * res[3] = A[N[6:7]]
 * res[4] = A[4 + N[0:1]]
 * res[5] = A[4 + N[2:3]]
 * res[6] = A[4 + N[4:5]]
 * res[7] = A[4 + N[6:7]]
 *
 * \tparam N
 * \param A
 * \return __m256i
 */
template <int N>
static __m256i permute256x32(__m256i A) {
  return _mm256_shuffle_epi32(A, N);
}
/**
 * \brief Permute given vector. The operation is equivalent to:
 * \brief res[0] = A[N[0:1]]
 * \brief res[1] = A[N[2:3]]
 *  res[2] = A[N[4:5]]
 *  res[3] = A[N[6:7]]
 *
 * \tparam N
 * \param A
 * \return __m128i
 */
template <int N>
static __m128i permute128x32(__m128i A) {
  return _mm_shuffle_epi32(A, N);
}
static __m256i mulmod256x32(__m256i A, __m256i B, __m256i M, __m256i M_Inv) {
  const __m256i prod0246 = _mm256_mul_epu32(A, B);
  const __m256i prod1357 =
      _mm256_mul_epu32(_mm256_srli_si256(A, 4), _mm256_srli_si256(B, 4));
  const __m256i qlo = _mm256_mul_epu32(prod0246, M_Inv);
  const __m256i qhi = _mm256_mul_epu32(prod1357, M_Inv);
  const __m256i tmplo = _mm256_mul_epu32(qlo, M);
  const __m256i tmphi = _mm256_mul_epu32(qhi, M);
  const __m256i alo = _mm256_sub_epi32(prod0246, tmplo);
  const __m256i ahi = _mm256_sub_epi32(prod1357, tmphi);
  const __m256i tmpd = shuffle256x32<_MM_SHUFFLE(3, 1, 3, 1)>(alo, ahi);
  const __m256i d = _mm256_shuffle_epi32(tmpd, _MM_SHUFFLE(3, 1, 2, 0));
  return _mm256_add_epi32(d, _mm256_and_si256(_mm256_srai_epi32(d, 31), M));
}
static __m256i reduce256x32(__m256i A, __m256i M, __m256i M_Inv) {
  const __m256i qlo = _mm256_mul_epu32(A, M_Inv);
  const __m256i qhi = _mm256_mul_epu32(_mm256_srli_si256(A, 4), M_Inv);
  const __m256i prodlo = _mm256_mul_epu32(qlo, M);
  const __m256i prodhi = _mm256_mul_epu32(qhi, M);
  const __m256i tmplo = _mm256_sub_epi32(_mm256_setzero_si256(), prodlo);
  const __m256i tmphi = _mm256_sub_epi32(_mm256_setzero_si256(), prodhi);
  const __m256i tmpd = shuffle256x32<_MM_SHUFFLE(3, 1, 3, 1)>(tmplo, tmphi);
  const __m256i d = _mm256_shuffle_epi32(tmpd, _MM_SHUFFLE(3, 1, 2, 0));
  return _mm256_add_epi32(d, _mm256_and_si256(_mm256_srai_epi32(d, 31), M));
}
static __m256i addmod256x32(__m256i A, __m256i B, __m256i M) {
  const __m256i d = _mm256_sub_epi32(_mm256_add_epi32(A, B), M);
  return _mm256_add_epi32(d, _mm256_and_si256(_mm256_srai_epi32(d, 31), M));
}
static __m256i submod256x32(__m256i A, __m256i B, __m256i M) {
  const __m256i d = _mm256_sub_epi32(A, B);
  return _mm256_add_epi32(d, _mm256_and_si256(_mm256_srai_epi32(d, 31), M));
}
namespace simd {
/**
 * \brief Mix and shuffle two vectors. The operation is equivalent to:
 * res[0] = A[0 + N[0:1]]
 * res[1] = A[0 + N[2:3]]
 * res[2] = B[0 + N[4:5]]
 * res[3] = B[0 + N[6:7]]
 * res[4] = A[4 + N[0:1]]
 * res[5] = A[4 + N[2:3]]
 * res[6] = B[4 + N[4:5]]
 * res[7] = B[4 + N[6:7]]
 *
 * \tparam N
 * \param A
 * \param B
 * \return __m256i
 */
template <int N>
static __m256i shuffle_x32(__m256i A, __m256i B) {
  return _mm256_castps_si256(
      _mm256_shuffle_ps(_mm256_castsi256_ps(A), _mm256_castsi256_ps(B), N));
}
/**
 * \brief Shuffle given two vectors (64 bit base). The operation is equivalent
 * to:
 *
 * res[0] = A[N[0]]
 * res[1] = B[N[1]]
 * res[2] = A[2 + N[2]]
 * res[3] = B[2 + N[3]]
 *
 * \tparam N
 * \param A
 * \param B
 * \return __m256i
 */
template <int N>
static __m256i shuffle_x64(__m256i A, __m256i B) {
  return _mm256_castpd_si256(
      _mm256_shuffle_pd(_mm256_castsi256_pd(A), _mm256_castsi256_pd(B), N));
}
/**
 * \brief Mix and shuffle given two vectors. The operation is equivalent to:
 *  Let D = A + B (concatenation of A and B).
 *  res[0] = D[N[0:1]]
 *  res[1] = D[N[2:3]]
 *
 * \tparam N
 * \param A
 * \param B
 * \return __m256i
 */
template <int N>
static __m256i shuffle_x128(__m256i A, __m256i B) {
  constexpr int Lo = (N >> 0) & 0b11;
  constexpr int Hi = (N >> 2) & 0b11;
  return _mm256_permute2f128_si256(A, B, (Hi << 4) | Lo);
}
/**
 * \brief Permute given vector. The operation is equivalent to:
 * res[0] = A[N[0:1]]
 * res[1] = A[N[2:3]]
 * res[2] = A[N[4:5]]
 * res[3] = A[N[6:7]]
 * res[4] = A[4 + N[0:1]]
 * res[5] = A[4 + N[2:3]]
 * res[6] = A[4 + N[4:5]]
 * res[7] = A[4 + N[6:7]]
 *
 * \tparam N
 * \param A
 * \return __m256i
 */
template <int N>
static __m256i permute_x32(__m256i A) {
  return _mm256_shuffle_epi32(A, N);
}
/**
 * \brief Permute given vector(64-bit base). The operation is equivalent to:
 * res[0] = A[N[0:1]]
 * res[1] = A[N[2:3]]
 * res[2] = A[N[4:5]]
 * res[3] = A[N[6:7]]
 *
 * \tparam N
 * \param A
 * \return __m256i
 */
template <int N>
static __m256i permute_x64(__m256i A) {
  return _mm256_permute4x64_epi64(A, N);
}
/**
 * \brief Permute given vector(128-bit base). The operation is equivalent to:
 *  res[0] = A[N[0]]
 *  res[1] = A[N[1]]
 *
 * \tparam N
 * \param A
 * \return __m256i
 */
template <int N>
static __m256i permute_x128(__m256i A) {
  constexpr int Lo = ((N >> 0) & 1) ? 0b1110 : 0b0100;
  constexpr int Hi = ((N >> 1) & 1) ? 0b1110 : 0b0100;
  return _mm256_permute4x64_epi64(A, (Hi << 4) | Lo);
}
/**
 * \brief Perform A[I] = f (32-bit base)
 *
 * \tparam I
 * \param A
 * \param f
 * \return __m256i
 */
template <int I>
static __m256i insert_x32(__m256i A, uint32_t f) {
  return _mm256_insert_epi32(A, f, I);
}
/**
 * \brief Perform A[I] = f (64-bit base)
 *
 * \tparam I
 * \param A
 * \param f
 * \return __m256i
 */
template <int I>
static __m256i insert_x64(__m256i A, uint64_t f) {
  return _mm256_insert_epi64(A, f, I);
}
/**
 * \brief Extract ith element of A (32-bit base)
 *
 * \tparam I
 * \param A
 * \return uint32_t
 */
template <int I>
static uint32_t extract_x32(__m256i A) {
  return static_cast<uint32_t>(_mm256_extract_epi32(A, I));
}
/**
 * \brief Extract ith element of A (64-bit base)
 *
 * \tparam I
 * \param A
 * \return uint64_t
 */
template <int I>
static uint64_t extract_x64(__m256i A) {
  return static_cast<uint64_t>(_mm256_extract_epi64(A, I));
}
/**
 * \brief For each i, perform A[i] = N[i] ? -A[i] : A[i];
 *
 * \tparam N
 * \param A
 * \return __m256i
 */
template <int N>
static __m256i mask_negate(__m256i A) {
  const auto mask = _mm256_set_epi32(-((N >> 7) & 1) | 1, -((N >> 6) & 1) | 1,
                                     -((N >> 5) & 1) | 1, -((N >> 4) & 1) | 1,
                                     -((N >> 3) & 1) | 1, -((N >> 2) & 1) | 1,
                                     -((N >> 1) & 1) | 1, -((N >> 0) & 1) | 1);
  return _mm256_sign_epi32(A, mask);
}
}  // namespace simd
}  // namespace algo::detail