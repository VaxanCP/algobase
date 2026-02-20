#pragma once

#include <immintrin.h>

#include <bit>
#include <limits>

#include "./def.hpp"
#include "./typing.hpp"
/*
@internal/base/typing.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo::detail {
/**
 * \brief Return the number of leading zeros in n.
 *
 * \tparam Tp integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int count_lz(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return std::countl_zero<Up>(n);
}
/**
 * \brief Return the number of trailing zeros in n
 *
 * \tparam Tp integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int count_tz(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return std::countr_zero<Up>(n);
}
/**
 * \brief Return the number of bits needed to represent n. If n is zero,
 * return zero
 *
 * \tparam Tp integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int bit_width(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return std::bit_width<Up>(n);
}
/**
 * \brief Return largest 2^k not greater than n
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp floor_pow2(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return static_cast<Tp>(std::bit_floor<Up>(n));
}
/**
 * \brief return smallest 2^k not less than n
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp ceil_pow2(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return static_cast<Tp>(std::bit_ceil<Up>(n));
}
/**
 * \brief Return [log2(n)]. If n is zero, return -1
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int floor_log2(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  constexpr int Nd = std::numeric_limits<Up>::digits;
  return Nd - 1 - std::countl_zero<Up>(n);
}
/**
 * \brief Return ceil(log2(n))
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int ceil_log2(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  constexpr int Nd = std::numeric_limits<Up>::digits;
  if (n == 0) { return 0; }
  return Nd - std::countl_zero<Up>(n - 1);
}
/**
 * \brief Copy the bits in n to result except for lowest set bits in n.
 * Equivalently, return n & (n - 1)
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp blsr(Tp n) {
  if (std::is_constant_evaluated()) { return n & (n - 1); }
  if constexpr (dword_fittable<Tp>) {
    return static_cast<Tp>(_blsr_u32(n));
  } else if constexpr (qword_fittable<Tp>) {
    return static_cast<Tp>(_blsr_u64(n));
  } else {
    return n & (n - 1);
  }
}
/**
 * \brief Extract the lowest set bit from n. Equivalently, return n & -n
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <typename Tp>
constexpr Tp blsi(Tp n) {
  if (std::is_constant_evaluated()) { return n & -n; }
  if constexpr (dword_fittable<Tp>) {
    return static_cast<Tp>(_blsi_u32(n));
  } else if constexpr (qword_fittable<Tp>) {
    return static_cast<Tp>(_blsi_u64(n));
  } else {
    return n & -n;
  }
}
/**
 * \brief Return the number of set bits in n
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr int
 */
template <typename Tp>
constexpr int popcount(Tp n) {
  using Up = std::make_unsigned_t<Tp>;
  return std::popcount<Up>(n);
}
}  // namespace algo::detail