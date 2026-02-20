#pragma once

#include <array>
#include <cstdint>
#include <span>

#include "./base/bit-base.hpp"
#include "./base/typing.hpp"
#include "./modular.hpp"

/*
@internal/base/bit-base.hpp
@internal/base/typing.hpp
@internal/modular.hpp
*/

// makecode
namespace algo::detail {
// Scaler version of discrete fourier transform implementations.
template <typename ModNum>
struct dif_info {
  static constexpr uint32_t Mod = ModNum::mod();      // Modulus
  static constexpr uint32_t Cap = count_tz(Mod - 1);  //
  static constexpr uint32_t Cof = (Mod - 1) >> Cap;
  static constexpr uint32_t Gen = primitive_root_prime(Mod);
  std::array<ModNum, Cap + 1> root;
  constexpr dif_info() noexcept {
    // g^(c)
    auto rot = ModNum(Gen).pow(Cof);
    // root[i] = g^(c*(2^(k-i))) (Which is root of unity of order 2^i)
    for (int i = Cap; i >= 0; --i) {
      // In each step, rot = root[i], so rot*rot = g^(c*(2^(k-i+1))) = root[i-1]
      root[i] = rot;
      rot *= rot;
    }
  }
};
template <typename ModNum>
struct dit_info {
  static constexpr uint32_t Mod = ModNum::mod();      // Modulus
  static constexpr uint32_t Cap = count_tz(Mod - 1);  //
  static constexpr uint32_t Cof = (Mod - 1) >> Cap;
  static constexpr uint32_t Gen = primitive_root_prime(Mod);
  std::array<ModNum, Cap + 1> root;
  constexpr dit_info() noexcept {
    // 1/g^(c)
    auto rot = ModNum(Gen).pow(Cof).inv();
    // root[i] = g^(c*(2^(k-i))) (Which is root of unity of order 2^i)
    for (int i = Cap; i >= 0; --i) {
      // In each step, rot = root[i], so rot*rot = g^(c*(2^(k-i+1))) = root[i-1]
      root[i] = rot;
      rot *= rot;
    }
  }
};
/**
 * \brief Do inverse fourier transform and store the result in-place. out.size()
 * must be a power of 2.
 *
 * \tparam ModNum
 * \param out
 */
template <typename ModNum>
void dif_butterfly(std::span<ModNum> out) {
  static const dif_info<ModNum> info{};
  // assert
  static_assert(info.Mod % 4 == 1, "Bad modulus for NTT");
  // Imaginary unit
  const ModNum imag = info.root[2];
  // The number of points
  const int n = static_cast<int>(out.size());
  // The number of stages
  const int k = floor_log2(n);
  // Hint
  assume(n > 0);
  // Do radix-4 butterfly until the number of points is not greater than n.
  for (int h = 0; h + (k % 2) < k; h += 2) {
    // 2^(k-h) is the current number of points. (Which is greater or equal to 4)
    const int m = 1 << (k - h);
    // Hint
    assume(m > 0);
    assume(m % 4 == 0);
    // w^1
    const ModNum wn_2 = info.root[k - h];
    // w^2
    const ModNum wn_1 = wn_2 * wn_2;
    // w^3
    const ModNum wn_3 = wn_1 * wn_2;
    for (int j = 0; j < n; j += m) {
      // Set loop invariants
      ModNum w1 = 1;
      ModNum w2 = 1;
      ModNum w3 = 1;
      // Do radix-4 butterfly
      for (int i = 0; i < m / 4; ++i) {
        /*
        out[i + j + m / 4 * 0] = (x0 +     x1 + x2 +     x3) * w^0;
        out[i + j + m / 4 * 1] = (x0 -     x1 + x2 -     x3) * w^2;
        out[i + j + m / 4 * 2] = (x0 + j * x1 - x2 - j * x3) * w^1;
        out[i + j + m / 4 * 3] = (x0 - j * x1 - x2 + j * x3) * w^3;
        */
        // Refer to: https://hackmd.io/@akshayk07/ryn-yR7qr
        const ModNum x0 = out[i + j + m / 4 * 0];
        const ModNum x1 = out[i + j + m / 4 * 1];
        const ModNum x2 = out[i + j + m / 4 * 2];
        const ModNum x3 = out[i + j + m / 4 * 3];
        /*
        (x0 + x2) + (x1 + x3)
        (x0 + x2) - (x1 + x3)
        (x0 - x2) + j*(x1 - x3)
        (x0 - x2) - j*(x1 - x3)
        */
        // x0 + x2
        const ModNum x0plx2 = x0 + x2;
        // x0 - x2
        const ModNum x0nex2 = x0 - x2;
        // x1 + x3
        const ModNum x1plx3 = x1 + x3;
        // j*(x1 - x3)
        const ModNum x1nex3j = (x1 - x3) * imag;
        // x0 + x1 + x2 + x3
        const ModNum e0 = x0plx2 + x1plx3;
        // x0 + x2 - x1 - x3
        const ModNum e1 = x0plx2 - x1plx3;
        // x0 - x2 + j*(x1-x3)
        const ModNum e2 = x0nex2 + x1nex3j;
        // x0-x2-j*(x1-x3)
        const ModNum e3 = x0nex2 - x1nex3j;
        // Store the result
        // (x0 + x2 + x1 + x3) * 1
        out[i + j + m / 4 * 0] = e0;
        // (x0 + x2 - x1 - x3) * w1
        out[i + j + m / 4 * 1] = e1 * w1;
        // (x0 - x2 + j*(x1 - x3))*w2
        out[i + j + m / 4 * 2] = e2 * w2;
        // (x0 - x2 - j*(x1-x3))*w3
        out[i + j + m / 4 * 3] = e3 * w3;
        // Update loop invariants
        w1 *= wn_1;
        w2 *= wn_2;
        w3 *= wn_3;
      }
    }
  }
  // When k is odd, do one radix-2 butterfly
  if (k % 2 == 1) {
    for (int i = 0; i < n; i += 2) {
      /*
        const mint a0 = out[i + j];
        const mint a1 = out[i + j + m / 2];
        out[i + j] = a0 + a1;
        out[i + j + m / 2] = (a0 - a1) * w;
      */
      const ModNum x0 = out[i + 0];
      const ModNum x1 = out[i + 1];
      out[i + 0] = x0 + x1;
      out[i + 1] = x0 - x1;
    }
  }
}
/**
 * \brief Do inverse fourier transform using decimation in time method then
 * store the result in-place.
 *
 * \tparam ModNum
 * \param out
 */
template <typename ModNum>
void dit_butterfly(std::span<ModNum> out) {
  // Setup
  static const dit_info<ModNum> info{};
  // assert
  static_assert(info.Mod % 4 == 1, "Bad modulus for NTT");
  // Imaginary unit
  const ModNum imag = info.root[2];
  // The number of points
  const int n = static_cast<int>(out.size());
  // The number of stages
  const int k = floor_log2(n);
  // If k is odd, do one pass radix-2 butterfly to resolve oddity issue.
  if (k % 2 == 1) {
    for (int i = 0; i < n; i += 2) {
      const ModNum x0 = out[i + 0];
      const ModNum x1 = out[i + 1];
      out[i + 0] = x0 + x1;
      out[i + 1] = x0 - x1;
    }
  }
  // Do radix-4 to cover remaining stages
  for (int h = 2 + (k % 2); h <= k; h += 2) {
    // 2^h is the current number of stages
    const int m = 1 << h;
    // Hint for the compiler
    assume(m > 0);
    assume(m % 4 == 0);
    // w^1
    const ModNum wn_2 = info.root[h];
    // w^2
    const ModNum wn_1 = wn_2 * wn_2;
    // w^3
    const ModNum wn_3 = wn_1 * wn_2;
    for (int j = 0; j < n; j += m) {
      // Loop variables
      ModNum w1 = 1;
      ModNum w2 = 1;
      ModNum w3 = 1;
      for (int i = 0; i < m / 4; ++i) {
        /*
          const mint x0 = out[base + m / 4 * 0] * w^0;
          const mint x1 = out[base + m / 4 * 1] * w^2;
          const mint x2 = out[base + m / 4 * 2] * w^1;
          const mint x3 = out[base + m / 4 * 3] * w^3;
        */
        /*
          out[base + m / 4 * 0] = x0 + x1 +     x2 +     x3;
          out[base + m / 4 * 1] = x0 - x1 + j * x2 - j * x3;
          out[base + m / 4 * 2] = x0 + x1 -     x2 -     x3;
          out[base + m / 4 * 3] = x0 - x1 - j * x2 + j * x3;
         */
        // Refer to:http://dsp-book.narod.ru/FFTBB/0270_PDF_C11.pdf
        const ModNum x0 = out[i + j + m / 4 * 0];
        const ModNum x1 = out[i + j + m / 4 * 1] * w1;
        const ModNum x2 = out[i + j + m / 4 * 2] * w2;
        const ModNum x3 = out[i + j + m / 4 * 3] * w3;
        // x0 + x1
        const ModNum x0plx1 = x0 + x1;
        // x0 - x1
        const ModNum x0nex1 = x0 - x1;
        // x2 + x3
        const ModNum x2plx3 = x2 + x3;
        // x2 - x3
        const ModNum x2nex3j = (x2 - x3) * imag;
        // Store the results
        // (x0 + x1) + (x2 + x3)
        out[i + j + m / 4 * 0] = x0plx1 + x2plx3;
        // (x0 - x1) + j*(x2 - x3)
        out[i + j + m / 4 * 1] = x0nex1 + x2nex3j;
        // (x0 + x1) - (x2 + x3)
        out[i + j + m / 4 * 2] = x0plx1 - x2plx3;
        // (x0 - x1) - j*(x2 - x3)
        out[i + j + m / 4 * 3] = x0nex1 - x2nex3j;
        // Update loop invariants
        w1 *= wn_1;
        w2 *= wn_2;
        w3 *= wn_3;
      }
    }
  }
}
}  // namespace algo::detail