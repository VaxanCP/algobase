#pragma once

#include <immintrin.h>

#include <cassert>
#include <numeric>
#include <span>

#include "./helpers.hpp"
#include "./modint.hpp"

/*
@internal/simd/helpers.hpp
@internal/simd/modint.hpp
*/

// makecode
namespace algo::detail::simd {
/*
  This namespace is dedicated to precompute constants needed for fast
  fourier transform implementation. Let M be a prime number with M =
  c*2^k+1 and w0 = g^c where g is a primitive root of M (clealy, g exists).
  Then if we set w[x]
  := g^(c*2^(k-x)) (0 <= x <= k),it follows that 1.(w[x])^(2^k) =
  g^(c*2^(k)) = 1 mod M 2.(w[x])^m != 1 mod M (for m != 2^x) Therefore,
  w[x] is root of unity of 2^x under modulo M. In addition, w[x]^2 is also
  a root of unity of 2^(x-1). Thus, for n = 2^l (where 0<=l<=k), we can
  compute fast fourier transform in O(n*Logn).
*/

/**
 * \brief Information needed for decimation in frequency DFT.
 *
 * \tparam MOD
 */
template <uint32_t MOD>
struct dif_info {
  using mint_type = montgomery_modint<MOD>;
  using mintx8_type = montgomery_modintx8<MOD>;
  static constexpr uint32_t Cap = count_tz(MOD - 1);
  static constexpr uint32_t Cof = (MOD - 1) >> Cap;
  static constexpr uint32_t Gen = primitive_root_prime(MOD);
  // root[i]^(2^i) = 1
  std::array<mintx8_type, Cap + 1> rootx8;
  std::array<mintx8_type, Cap + 1> stepx8;
  dif_info() noexcept {
    // root[i] = g^(c*(2^(k-i))) then root[i]^(2^i)=1
    auto rot = mint_type(Gen).pow(Cof);
    // rot^(2^i)=1 for each step i.
    for (int i = Cap + 1; i-- > 0;) {
      mint_type f = 1;
      // buf = [rot^0 rot^1 rot^2 rot^3 ... rot^7]
      alignas(__m256i) std::array<mint_type, 8> buf;
      for (int j = 0; j < 8; ++j) {
        // For each step j, f = rot^j.
        buf[j] = f;
        f *= rot;
      }
      // [rot^0 rot^1 rot^2...]
      rootx8[i].load_aligned(buf.data());
      // f = rot^8
      buf.fill(f);
      // [w^8 w^8 w^8...]
      stepx8[i].load_aligned(buf.data());
      // Update loop invariants
      rot *= rot;
    }
  }
};
/**
 * \brief Information metaclass needed for decimation in time DFT.
 *
 * \tparam MOD
 */
template <uint32_t MOD>
struct dit_info {
  using mint_type = montgomery_modint<MOD>;
  using mintx8_type = montgomery_modintx8<MOD>;
  static constexpr uint32_t Cap = count_tz(MOD - 1);
  static constexpr uint32_t Cof = (MOD - 1) >> Cap;
  static constexpr uint32_t Gen = primitive_root_prime(MOD);
  std::array<mintx8_type, Cap + 1> irootx8;
  std::array<mintx8_type, Cap + 1> istepx8;
  dit_info() noexcept {
    // root[i] = g^(c*(2^(k-i))) then root[i]^(2^i)=1
    // Note that rot is not in the montgomery space.
    auto rot = mint_type(Gen).pow(Cof).inv();
    // rot^(2^i)=1 for each step i.
    for (int i = Cap + 1; i-- > 0;) {
      mint_type f = 1;
      // buf = [rot^0 rot^1 rot^2 rot^3 ... rot^7]
      alignas(__m256i) std::array<mint_type, 8> buf;
      for (int j = 0; j < 8; ++j) {
        // For each step j, f = rot^j.
        buf[j] = f;
        f *= rot;
      }
      // [rot^0 rot^1 rot^2...]
      irootx8[i].load_aligned(buf.data());
      // f = rot^8
      buf.fill(f);
      // [w^8 w^8 w^8...]
      istepx8[i].load_aligned(buf.data());
      // Update loop invariants
      rot *= rot;
    }
  }
};
/**
 * \brief Do discrete fourier transform using decimation in frequency method and
 * store the result in-place. The size of input must be a power of 2.
 *
 * \tparam MOD
 * \param out
 */
template <uint32_t MOD = 998244353>
void dif_butterfly(std::span<montgomery_modintx8<MOD>> out) {
  static_assert(MOD % 8 == 1, "Bad modulus for NTT");
  using mint_type = montgomery_modint<MOD>;
  using mintx8_type = montgomery_modintx8<MOD>;
  // Set up
  static const dif_info<MOD> info{};
  // root[3]^0 ( = j^(0))
  const uint32_t rot3_0 =
      mint_type(info.Gen).pow(0 * (info.Cof << (info.Cap - 3))).get();
  // root[3]^1 ( = j^(1/2))
  const uint32_t rot3_1 =
      mint_type(info.Gen).pow(1 * (info.Cof << (info.Cap - 3))).get();
  // root[3]^2 ( = j^(1))
  const uint32_t rot3_2 =
      mint_type(info.Gen).pow(2 * (info.Cof << (info.Cap - 3))).get();
  // root[3]^3 ( = j^(3/2))
  const uint32_t rot3_3 =
      mint_type(info.Gen).pow(3 * (info.Cof << (info.Cap - 3))).get();
  // Register setup
  const mintx8_type imagx8{_mm256_set1_epi32(rot3_2)};
  // The size of inputs divided by 8
  const int n = static_cast<int>(out.size());
  // The number of stages
  const int k = floor_log2(n) + 3;
  // n > 0
  assume(n > 0);
  // First do radix-4 butterfly while the number of points greater than 16 to
  // utilize the vectorization, then compute remaining stages afterward.
  for (int h = 0; h + 4 < k; h += 2) {
    // 2^(k-h) is the current number of points. Since k-h>4, the number of
    // points >= 32. The current number of points divided by 8
    const int m = 1 << (k - h - 3);
    // Hint
    assume(m > 0);
    assume(m % 4 == 0);
    // w denotes the root of unity of order m
    // [w^0,w^1,w^2,w^3,w^4,w^5,w^6,w^7]
    const auto wn2x8 = info.rootx8[k - h];
    // [w^0,w^2,w^4,w^6,w^8,w^10,w^12,w^14]
    const auto wn1x8 = wn2x8 * wn2x8;
    // [w^0,w^3,w^6,w^9,w^12,w^15,w^18,w^21]
    const auto wn3x8 = wn1x8 * wn2x8;
    // The stride to keep the loop invariants
    // [w^8,w^8,...]
    const auto sw2x8 = info.stepx8[k - h];
    // [w^16,w^16,..]
    const auto sw1x8 = sw2x8 * sw2x8;
    // [w^24,w^24,...]
    const auto sw3x8 = sw1x8 * sw2x8;
    for (int j = 0; j < n; j += m) {
      // Set loop invariants
      auto w1x8 = wn1x8;
      auto w2x8 = wn2x8;
      auto w3x8 = wn3x8;
      // Do radix-4 butterfly
      for (int i = 0; i < m / 4; ++i) {
        /*
        out[i + j + m / 4 * 0] = (x0 +     x1 + x2 +     x3) * w^0;
        out[i + j + m / 4 * 1] = (x0 -     x1 + x2 -     x3) * w^2;
        out[i + j + m / 4 * 2] = (x0 + j * x1 - x2 - j * x3) * w^1;
        out[i + j + m / 4 * 3] = (x0 - j * x1 - x2 + j * x3) * w^3;
        */
        // Refer to: https://hackmd.io/@akshayk07/ryn-yR7qr
        // We need to divide indexes by 8 because data type is packed 8
        // integers.
        const auto x0 = out[i + j + m / 4 * 0];
        const auto x1 = out[i + j + m / 4 * 1];
        const auto x2 = out[i + j + m / 4 * 2];
        const auto x3 = out[i + j + m / 4 * 3];
        /*
        (x0 + x2) + (x1 + x3)
        (x0 + x2) - (x1 + x3)
        (x0 - x2) + j*(x1 - x3)
        (x0 - x2) - j*(x1 - x3)
        */
        // x0 + x2
        const auto x0plx2 = x0 + x2;
        // x0 - x2
        const auto x0nex2 = x0 - x2;
        // x1 + x3
        const auto x1plx3 = x1 + x3;
        // j*(x1 - x3)
        const auto x1nex3j = (x1 - x3) * imagx8;
        // x0 + x1 + x2 + x3
        const auto e0 = x0plx2 + x1plx3;
        // x0 + x2 - x1 - x3
        const auto e1 = x0plx2 - x1plx3;
        // x0 - x2 + j*(x1-x3)
        const auto e2 = x0nex2 + x1nex3j;
        // x0-x2-j*(x1-x3)
        const auto e3 = x0nex2 - x1nex3j;
        // Store the result
        // (x0 + x2 + x1 + x3) * 1
        out[i + j + m / 4 * 0] = e0;
        // (x0 + x2 - x1 - x3) * w1
        out[i + j + m / 4 * 1] = e1 * w1x8;
        // (x0 - x2 + j*(x1 - x3))*w2
        out[i + j + m / 4 * 2] = e2 * w2x8;
        // (x0 - x2 - j*(x1-x3))*w3
        out[i + j + m / 4 * 3] = e3 * w3x8;
        // Update loop invariants
        w1x8 *= sw1x8;
        w2x8 *= sw2x8;
        w3x8 *= sw3x8;
      }
    }
  }
  // Remaining stages are either 3 or 4 depending on the oddity of k.
  if (k % 2 == 0) {
    // Remaining stages are 4, so do one radix-2 butterfly then do radix-3
    // butterfly to finish.
    /*
      const mint a0 = out[i + j];
      const mint a1 = out[i + j + m / 2];
      out[i + j] = a0 + a1;
      out[i + j + m / 2] = (a0 - a1) * w;
    */
    // m = 16 = 2^4
    const auto wnx8 = info.rootx8[4];
    for (int j = 0; j < n; j += 2) {
      // Divide index by 8
      const auto x0 = out[j + 0];
      const auto x1 = out[j + 1];
      // x0 + x1
      out[j + 0] = x0 + x1;
      // w*(x0-x1)
      out[j + 1] = (x0 - x1) * wnx8;
    }
  }
  // Last radix-8
  {
    // [1 1 1 1 1 1 j j]
    const auto v2x8 = _mm256_set_epi32(rot3_2, rot3_2, rot3_0, rot3_0, rot3_0,
                                       rot3_0, rot3_0, rot3_0);
    // [1 1 1 j 1 j^(1/2) 1 j^(3/2)]
    const auto v3x8 = _mm256_set_epi32(rot3_3, rot3_0, rot3_1, rot3_0, rot3_2,
                                       rot3_0, rot3_0, rot3_0);
    for (int j = 0; j < n; ++j) {
      // [x0 x1 x2 x3 x4 x5 x6 x7]
      const auto x1 = out[j].ymm;
      // [x4 x5 x6 x7 x0 x1 x2 x3]
      const auto y1 = permute_x64<_MM_SHUFFLE(1, 0, 3, 2)>(x1);
      // Negate last half ([x0 x1 x2 x3 -x4 -x5 -x6 -x7])
      const auto z1 = mask_negate<0b11110000>(x1);
      // Normalize z1 [x0 x1 x2 x3 -x4 -x5 -x6 -x7]
      const auto t1 = mintx8_type::normalize(z1);
      // [x0+x4 x1+x5 x2+x6 x3+x7 x0-x4 x1-x5 x2-x6 x3-x7]
      const auto w1 = mintx8_type::add_mod(y1, t1);
      // [x0+x4 x1+x5 x2+x6 x3+x7 x0-x4 x1-x5 j(x2-x6) j(x3-x7)]
      const auto x2 = mintx8_type::mul_mod(w1, v2x8);
      // [x2+x6 x3+x7 x0+x4 x1+x5 j(x3-x7) j(x2-x6) x0-x4 x1-x5]
      const auto y2 = permute_x32<_MM_SHUFFLE(1, 0, 3, 2)>(x2);
      // [x0+x4 x1+x5 -(x2+x6) -(x3+x7) x0-x4 x1-x5 -j(x2-x6) -j(x3-x7)]
      const auto z2 = mask_negate<0b11001100>(x2);
      // Normalize z2
      const auto t2 = mintx8_type::normalize(z2);
      // [(x0+x4)+(x2+x6) (x1+x5)+(x3+x7) (x0+x4)-(x2+x6) (x1+x5)-(x3+x7)
      // (x0-x4)+j(x3-x7) (x1-x5)+j(x2-x6) (x0-x4)-j(x2-x6) (x1-x5)-j(x3-x7)]
      const auto w2 = mintx8_type::add_mod(y2, t2);
      // [(x0+x4)+(x2+x6) (x1+x5)+(x3+x7) (x0+x4)-(x2+x6) j((x1+x5)-(x3+x7))
      // (x0-x4)+j(x3-x7) j^(1/2)((x1-x5)+j(x2-x6)) (x0-x4)-j(x2-x6)
      // j^(3/2)((x1-x5)-j(x3-x7))]
      const auto x3 = mintx8_type::mul_mod(w2, v3x8);
      // [(x1+x5)+(x3+x7) (x0+x4)+(x2+x6) j((x1+x5)-(x3+x7)) (x0+x4)-(x2+x6)
      // j^(1/2)((x1-x5)+j(x2-x6)) (x0-x4)+j(x3-x7) j^(3/2)((x1-x5)-j(x3-x7))
      // (x0-x4)-j(x2-x6)]
      const auto y3 = permute_x32<_MM_SHUFFLE(2, 3, 0, 1)>(x3);
      // [(x0+x4)+(x2+x6) -(x1+x5)-(x3+x7) (x0+x4)-(x2+x6) -j((x1+x5)-(x3+x7))
      // (x0-x4)+j(x3-x7) -j^(1/2)((x1-x5)+j(x2-x6)) (x0-x4)-j(x2-x6)
      // -j^(3/2)((x1-x5)-j(x3-x7))]
      const auto z3 = mask_negate<0b10101010>(x3);
      // Normalize z3
      const auto t3 = mintx8_type::normalize(z3);
      const auto w3 = mintx8_type::add_mod(y3, t3);
      out[j].ymm = w3;
    }
  }
}
/**
 * \brief Do Inverse fourier transform using decimation in time method then
 * store the result in-place.
 *
 * \tparam MOD
 * \param out
 */
template <uint32_t MOD = 998244353>
void dit_butterfly(std::span<montgomery_modintx8<MOD>> out) {
  static_assert(MOD % 8 == 1, "Bad modulus for NTT");
  using mint_type = montgomery_modint<MOD>;
  using mintx8_type = montgomery_modintx8<MOD>;
  // Set up
  static const dit_info<MOD> info{};
  // root[3]^0 ( = j^(0))
  const uint32_t rot3_0 =
      mint_type(info.Gen).pow(0 * (info.Cof << (info.Cap - 3))).inv().get();
  // root[3]^1 ( = j^(1/2))
  const uint32_t rot3_1 =
      mint_type(info.Gen).pow(1 * (info.Cof << (info.Cap - 3))).inv().get();
  // root[3]^2 ( = j^(1))
  const uint32_t rot3_2 =
      mint_type(info.Gen).pow(2 * (info.Cof << (info.Cap - 3))).inv().get();
  // root[3]^3 ( = j^(3/2))
  const uint32_t rot3_3 =
      mint_type(info.Gen).pow(3 * (info.Cof << (info.Cap - 3))).inv().get();
  // Register set up
  const mintx8_type imagx8{_mm256_set1_epi32(rot3_2)};
  // The size of inputs divided by 8
  const int n = static_cast<int>(out.size());
  // The number of stages
  const int k = floor_log2(n) + 3;
  // Hint
  assume(n > 0);
  // First do radix-8 butterfly
  {
    /*
      target:
      (x0+x1)+(x2+x3)+(x4+x5)+(x6+x7)
      (x0-x1)+j(x2-x3)+j^(1/2)((x4-x5)+j(x6-x7))
      (x0+x1)-(x2+x3)+j((x4+x5)-(x6+x7))
      (x0-x1)-j(x2-x3)+j^(3/2)((x4-x5)-j(x6-x7))
      (x0+x1)+(x2+x3)-(x4+x5)-(x6+x7)
      (x0-x1)+j(x2-x3)-j^(1/2)((x4-x5)+j(x6-x7))
      (x0+x1)-(x2+x3)-j((x4+x5)-(x6+x7))
      (x0-x1)-j(x2-x3)+j^(3/2)((x4-x5)-j(x6-x7))
    */
    // [1 1 1 j 1 1 1 j]
    const auto v2x8 = _mm256_set_epi32(rot3_2, rot3_0, rot3_0, rot3_0, rot3_2,
                                       rot3_0, rot3_0, rot3_0);
    // [1 1 1 1 1 j^(1/2) j j^(3/2)]
    const auto v3x8 = _mm256_set_epi32(rot3_3, rot3_2, rot3_1, rot3_0, rot3_0,
                                       rot3_0, rot3_0, rot3_0);
    for (int j = 0; j < n; ++j) {
      // [x0 x1 x2 x3 x4 x5 x6 x7]
      const auto x1 = out[j].ymm;
      // [x1 x0 x3 x2 x5 x4 x7 x6]
      const auto y1 = permute_x32<_MM_SHUFFLE(2, 3, 0, 1)>(x1);
      // Negate x
      // [x0 -x1 x2 -x3 x4 -x5 x6 -x7]
      const auto z1 = mask_negate<0b10101010>(x1);
      // Normalize z1
      const auto t1 = mintx8_type::normalize(z1);
      // [x0+x1 x0-x1 x2+x3 x2-x3 x4+x5 x4-x5 x6+x7 x6-x7]
      const auto w1 = mintx8_type::add_mod(y1, t1);
      // [x0+x1 x0-x1 x2+x3 j(x2-x3) x4+x5 x4-x5 x6+x7 j(x6-x7)]
      const auto x2 = mintx8_type::mul_mod(w1, v2x8);
      // [x2+x3 j(x2-x3) x0+x1 x0-x1 x6+x7 j(x6-x7) x4+x5 x4-x5]
      const auto y2 = permute_x32<_MM_SHUFFLE(1, 0, 3, 2)>(x2);
      // Negate x2
      // [x0+x1 x0-x1 -(x2+x3) -j(x2-x3) x4+x5 x4-x5 -(x6+x7) -j(x6-x7)]
      const auto z2 = mask_negate<0b11001100>(x2);
      // Normalize z2
      const auto t2 = mintx8_type::normalize(z2);
      // [(x0+x1)+(x2+x3) (x0-x1)+j(x2-x3) (x0+x1)-(x2+x3) (x0-x1)-j(x2-x3)
      // (x4+x5)+(x6+x7) (x4-x5)+j(x6-x7) (x4+x5)-(x6+x7) (x4-x5)-j(x6-x7)]
      const auto w2 = mintx8_type::add_mod(y2, t2);
      // [(x0+x1)+(x2+x3) (x0-x1)+j(x2-x3) (x0+x1)-(x2+x3) (x0-x1)-j(x2-x3)
      // (x4+x5)+(x6+x7) j^(1/2)((x4-x5)+j(x6-x7)) j((x4+x5)-(x6+x7))
      // j^(3/2)((x4-x5)-j(x6-x7))]
      const auto x3 = mintx8_type::mul_mod(w2, v3x8);
      // [(x4+x5)+(x6+x7) j^(1/2)((x4-x5)+j(x6-x7)) j((x4+x5)-(x6+x7))
      // j^(3/2)((x4-x5)-j(x6-x7)) (x0+x1)+(x2+x3) (x0-x1)+j(x2-x3)
      // (x0+x1)-(x2+x3) (x0-x1)-j(x2-x3)]
      const auto y3 = permute_x64<_MM_SHUFFLE(1, 0, 3, 2)>(x3);
      // Negate x
      // [(x0+x1)+(x2+x3) (x0-x1)+j(x2-x3) (x0+x1)-(x2+x3) (x0-x1)-j(x2-x3)
      // -((x4+x5)+(x6+x7)) -j^(1/2)((x4-x5)+j(x6-x7)) -j((x4+x5)-(x6+x7))
      // -j^(3/2)((x4-x5)-j(x6-x7))]
      const auto z3 = mask_negate<0b11110000>(x3);
      // Normalize z3
      const auto t3 = mintx8_type::normalize(z3);
      //  [(x0+x1)+(x2+x3)+(x4+x5)+(x6+x7)
      // (x0-x1)+j(x2-x3)+j^(1/2)((x4-x5)+j(x6-x7))
      // (x0+x1)-(x2+x3)+j((x4+x5)-(x6+x7))
      // (x0-x1)-j(x2-x3)+j^(3/2)((x4-x5)-j(x6-x7))
      // (x0+x1)+(x2+x3)-(x4+x5)-(x6+x7)
      // (x0-x1)+j(x2-x3)-j^(1/2)((x4-x5)+j(x6-x7))
      // (x0+x1)-(x2+x3)-j((x4+x5)-(x6+x7))
      // (x0-x1)-j(x2-x3)+j^(3/2)((x4-x5)-j(x6-x7))]
      const auto w3 = mintx8_type::add_mod(y3, t3);
      out[j].ymm = w3;
    }
  }
  if (k % 2 == 0) {
    // If k even, do one radix-2 butterfly
    // 16=2^4
    const auto wnx8 = info.irootx8[4];
    for (int j = 0; j < n; j += 2) {
      const auto x0 = out[j + 0];
      const auto x1 = out[j + 1] * wnx8;
      // x0 + x1
      out[j + 0] = x0 + x1;
      // x0 - x1
      out[j + 1] = x0 - x1;
    }
  }
  for (int h = 6 - (k % 2); h <= k; h += 2) {
    // current number of points divided by 8
    const int m = 1 << (h - 3);
    // Hint for the compiler
    assume(m > 0);
    assume(m % 4 == 0);
    // AVX vectorized radix-4 butterfly
    // [w^0,w^1,w^2,w^3,w^4,w^5,w^6,w^7]
    const auto wn2x8 = info.irootx8[h];
    // [w^0,w^2,w^4,w^6,w^8,w^10,w^12,w^14]
    const auto wn1x8 = wn2x8 * wn2x8;
    // [w^0,w^3,w^6,w^9,w^12,w^15,w^18,w^21]
    const auto wn3x8 = wn1x8 * wn2x8;
    // [w^8,w^8,...]
    const auto sw2x8 = info.istepx8[h];
    // [w^16,w^16,...]
    const auto sw1x8 = sw2x8 * sw2x8;
    // [w^24,w^24,...]
    const auto sw3x8 = sw1x8 * sw2x8;
    for (int j = 0; j < n; j += m) {
      // Set loop invariants
      auto w1x8 = wn1x8;
      auto w2x8 = wn2x8;
      auto w3x8 = wn3x8;
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
        const auto x0 = out[i + j + m / 4 * 0];
        const auto x1 = out[i + j + m / 4 * 1] * w1x8;
        const auto x2 = out[i + j + m / 4 * 2] * w2x8;
        const auto x3 = out[i + j + m / 4 * 3] * w3x8;
        // x0 + x1
        const auto x0plx1 = x0 + x1;
        // x0 - x1
        const auto x0nex1 = x0 - x1;
        // x2 + x3
        const auto x2plx3 = x2 + x3;
        // x2 - x3
        const auto x2nex3j = (x2 - x3) * imagx8;
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
        w1x8 *= sw1x8;
        w2x8 *= sw2x8;
        w3x8 *= sw3x8;
      }
    }
  }
}
}  // namespace algo::detail::simd