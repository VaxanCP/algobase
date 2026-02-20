#pragma once

#include "internal/base/bit-base.hpp"
/*
@internal/base/bit-base.hpp
*/
// makecode
namespace algo {
/**
 * \brief Return the number of leading zeros in n
 *
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int count_lz(Tp n) {
  return detail::count_lz(n);
}
/**
 * \brief Return the number of trailing zeros in n
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int count_tz(Tp n) {
  return detail::count_tz(n);
}
/**
 * \brief Return the number of bits needed to represent n
 *
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int bit_width(Tp n) {
  return detail::bit_width(n);
}
/**
 * \brief Return the largest power of 2 not greater than n
 *
 * \tparam Tp
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp floor_pow2(Tp n) {
  return detail::floor_pow2(n);
}
/**
 * \brief Return the smallest power of 2 not less than n
 *
 * \tparam Tp
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp ceil_pow2(Tp n) {
  return detail::ceil_pow2(n);
}
/**
 * \brief Return floor value of log2(n)
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int floor_log2(Tp n) {
  return detail::floor_log2(n);
}
/**
 * \brief Return ceil value of log2(n)
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int ceil_log2(Tp n) {
  return detail::ceil_log2(n);
}
/**
 * \brief Copy the bits in n to result except for lowest set bits in n.
 * Equivalently, return n & (n - 1)
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp blsr(Tp n) {
  return detail::blsr(n);
}
/**
 * \brief Extract the lowest set bit from n. Equivalently, return n & -n
 *
 * \tparam Tp Integer required
 * \param n
 * \return constexpr Tp
 */
template <detail::integer Tp>
constexpr Tp blsi(Tp n) {
  return detail::blsi(n);
}
/**
 * \brief Return the number of set bits in n
 * \tparam Tp
 * \param n
 * \return constexpr int
 */
template <detail::integer Tp>
constexpr int popcount(Tp n) {
  return detail::popcount(n);
}
}  // namespace algo