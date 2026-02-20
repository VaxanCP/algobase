#pragma once

#include <immintrin.h>

#include <cstdint>

// makecode
namespace algo::detail::simd {
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
}  // namespace algo::detail::simd