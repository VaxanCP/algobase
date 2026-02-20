#pragma once
#include <immintrin.h>

#include <array>
#include <bit>
#include <cstring>
#include <memory>
#include <new>

#include "../class/montgomery.hpp"
#include "./base/bit-base.hpp"
#include "./base/numeric-base.hpp"
#include "./base/typing.hpp"
#include "./modular.hpp"
#include "./simd.hpp"

/*
@internal/base/numeric-base.hpp
@internal/base/typing.hpp
@internal/modular.hpp
@internal/base/bit-base.hpp
@class/montgomery.hpp
@internal/simd.hpp
*/
// makecode
namespace algo::detail::fft {
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
// return zero-initialized storage
// note: size n is not bytes, but the number of ints
uint32_t* get_storage(size_t n) noexcept {
  constexpr std::align_val_t Algn = std::align_val_t(alignof(__m256i));
  void* mem = ::operator new[](n * sizeof(uint32_t), Algn, std::nothrow);
  mem = std::assume_aligned<alignof(__m256i)>(mem);
  memset(mem, 0, n * sizeof(uint32_t));
  return static_cast<uint32_t*>(mem);
}
void free_storage(uint32_t* ptr) noexcept {
  constexpr std::align_val_t Algn = std::align_val_t(alignof(__m256i));
  ::operator delete[](ptr, Algn);
}
void radix_permute(uint32_t* out, int n) {
  for (int i = 0, j = 0; i < n; ++i) {
    if (j < i) { std::swap(out[i], out[j]); }
    j ^= n - (n >> (count_tz(i + 1) + 1));
  }
}
// do dft and store the result in out
// precondition: n is power of two and sufficiently large
// precondition: out is aligned at 32 bytes
// precondition: values in out are already transformed into montogomery
// space note: resulting values are the values under montogomery space
template <uint32_t MOD, bool Inverse>
void do_butterfly(uint32_t* out, int n) {
  using mint_type = montgomery_modint<MOD>;
  constexpr uint32_t K = count_tz(MOD - 1);
  constexpr uint32_t C = (MOD - 1) >> K;
  constexpr uint32_t G = primitive_root_prime(MOD);
  constexpr static std::array<std::array<mint_type, 17>, K> Wn = [] {
    std::array<std::array<mint_type, 17>, K> res{};
    mint_type wn = mint_type(G).pow(C);
    if (Inverse) { wn = wn.inv(); }
    res[K - 1][0] = 1;
    for (uint32_t i = 1; i < 17; ++i) {
      res[K - 1][i] = res[K - 1][i - 1] * wn;
    }
    for (uint32_t i = K - 1; i-- > 0;) {
      for (uint32_t j = 0; j < 17; ++j) {
        res[i][j] = res[i + 1][j] * res[i + 1][j];
      }
    }
    return res;
  }();
  const __m256i reg_md = _mm256_set1_epi32(MOD);
  const __m256i reg_nd = _mm256_set1_epi32(mint_type::Nd);
  const __m256i reg_w1 = _mm256_set1_epi32(Wn[0][0].get());
  [[maybe_unused]] const __m256i reg_rsq = _mm256_set1_epi32(mint_type::R2);
  const __m256i reg_w2 = _mm256_set_epi32(
      Wn[1][1].get(), Wn[1][1].get(), Wn[1][0].get(), Wn[1][0].get(),
      Wn[1][1].get(), Wn[1][1].get(), Wn[1][0].get(), Wn[1][0].get());
  const __m256i reg_w3 = _mm256_set_epi32(
      Wn[2][3].get(), Wn[2][3].get(), Wn[2][2].get(), Wn[2][2].get(),
      Wn[2][1].get(), Wn[2][1].get(), Wn[2][0].get(), Wn[2][0].get());
  const __m256i reg_w4 = _mm256_set_epi32(
      Wn[3][7].get(), Wn[3][6].get(), Wn[3][5].get(), Wn[3][4].get(),
      Wn[3][3].get(), Wn[3][2].get(), Wn[3][1].get(), Wn[3][0].get());
  radix_permute(out, n);
  for (int i = 0; i < n; i += 16) {
    __m256i v0 = _mm256_load_si256((__m256i*)&out[i]);
    __m256i v1 = _mm256_load_si256((__m256i*)&out[i + 8]);
    if constexpr (!Inverse) {
      v0 = mulmod256x32(v0, reg_rsq, reg_md, reg_nd);
      v1 = mulmod256x32(v1, reg_rsq, reg_md, reg_nd);
    }
    __m256i u0 = shuffle256x32<_MM_SHUFFLE(2, 0, 2, 0)>(v0, v1);
    __m256i u1 = shuffle256x32<_MM_SHUFFLE(3, 1, 3, 1)>(v0, v1);
    u1 = mulmod256x32(u1, reg_w1, reg_md, reg_nd);
    v0 = addmod256x32(u0, u1, reg_md);
    v1 = submod256x32(u0, u1, reg_md);
    u0 = shuffle256x32<_MM_SHUFFLE(2, 0, 2, 0)>(v0, v1);
    u1 = shuffle256x32<_MM_SHUFFLE(3, 1, 3, 1)>(v0, v1);
    u1 = mulmod256x32(u1, reg_w2, reg_md, reg_nd);
    v0 = addmod256x32(u0, u1, reg_md);
    v1 = submod256x32(u0, u1, reg_md);
    u0 = shuffle256x128<0x20>(v0, v1);
    u1 = shuffle256x128<0x31>(v0, v1);
    u1 = mulmod256x32(u1, reg_w3, reg_md, reg_nd);
    v0 = addmod256x32(u0, u1, reg_md);
    v1 = submod256x32(u0, u1, reg_md);
    const __m256i tmp0 = shuffle256x32<_MM_SHUFFLE(2, 0, 2, 0)>(v0, v1);
    const __m256i tmp1 = shuffle256x32<_MM_SHUFFLE(3, 1, 3, 1)>(v0, v1);
    u0 = shuffle256x64<_MM_SHUFFLE(3, 1, 2, 0)>(tmp0);
    u1 = shuffle256x64<_MM_SHUFFLE(3, 1, 2, 0)>(tmp1);
    u1 = mulmod256x32(u1, reg_w4, reg_md, reg_nd);
    v0 = addmod256x32(u0, u1, reg_md);
    v1 = submod256x32(u0, u1, reg_md);
    _mm256_store_si256((__m256i*)&out[i], v0);
    _mm256_store_si256((__m256i*)&out[i + 8], v1);
  }
  for (int lg = 4, l = 32; l <= n; ++lg, l <<= 1) {
    const __m256i reg_stride = _mm256_set1_epi32(Wn[lg][16].get());
    const __m256i reg_fst = _mm256_set_epi32(
        Wn[lg][7].get(), Wn[lg][6].get(), Wn[lg][5].get(), Wn[lg][4].get(),
        Wn[lg][3].get(), Wn[lg][2].get(), Wn[lg][1].get(), Wn[lg][0].get());
    const __m256i reg_sec = _mm256_set_epi32(
        Wn[lg][15].get(), Wn[lg][14].get(), Wn[lg][13].get(), Wn[lg][12].get(),
        Wn[lg][11].get(), Wn[lg][10].get(), Wn[lg][9].get(), Wn[lg][8].get());
    for (int i = 0; i < n; i += l) {
      __m256i reg_wf = reg_fst;
      __m256i reg_ws = reg_sec;
      for (int j = 0; j < (l >> 1); j += 16) {
        __m256i v0 = _mm256_load_si256((__m256i*)&out[i + j]);
        __m256i v1 = _mm256_load_si256((__m256i*)&out[i + j + 8]);
        __m256i v2 = _mm256_load_si256((__m256i*)&out[i + j + (l >> 1)]);
        __m256i v3 = _mm256_load_si256((__m256i*)&out[i + j + (l >> 1) + 8]);
        v2 = mulmod256x32(v2, reg_wf, reg_md, reg_nd);
        v3 = mulmod256x32(v3, reg_ws, reg_md, reg_nd);
        _mm256_store_si256((__m256i*)&out[i + j], addmod256x32(v0, v2, reg_md));
        _mm256_store_si256((__m256i*)&out[i + j + 8],
                           addmod256x32(v1, v3, reg_md));
        _mm256_store_si256((__m256i*)&out[i + j + (l >> 1)],
                           submod256x32(v0, v2, reg_md));
        _mm256_store_si256((__m256i*)&out[i + j + (l >> 1) + 8],
                           submod256x32(v1, v3, reg_md));
        reg_wf = mulmod256x32(reg_wf, reg_stride, reg_md, reg_nd);
        reg_ws = mulmod256x32(reg_ws, reg_stride, reg_md, reg_nd);
      }
    }
  }
  if constexpr (Inverse) {
    const auto ninv = mint_type(static_cast<uint32_t>(n)).inv();
    const __m256i reg_ninv = _mm256_set1_epi32(ninv.get());
    for (int i = 0; i < n; i += 16) {
      __m256i v0 = _mm256_load_si256((__m256i*)&out[i]);
      __m256i v1 = _mm256_load_si256((__m256i*)&out[i + 8]);
      v0 = mulmod256x32(v0, reg_ninv, reg_md, reg_nd);
      v1 = mulmod256x32(v1, reg_ninv, reg_md, reg_nd);
      _mm256_store_si256((__m256i*)&out[i], reduce256x32(v0, reg_md, reg_nd));
      _mm256_store_si256((__m256i*)&out[i + 8],
                         reduce256x32(v1, reg_md, reg_nd));
    }
  }
}
// convolute two vector c1 and c2 under modulo MOD and store the result in
// c1 precondition: n is power of 2 and sufficiently large precondition: c1
// and c2 are both the same size (n) precondition: c1 and c2 are both
// aligned at 32 bytes precondition: values of both c1 and c2 are already
// transformed into montogomery space note: resulting values of c1 are
// values under montogomery space note: resulting values of c2 are
// undefined
template <uint32_t MOD>
void convolution(uint32_t* __restrict__ c1, uint32_t* __restrict__ c2, int n) {
  using mint_type = montgomery_modint<MOD>;
  const __m256i reg_md = _mm256_set1_epi32(MOD);
  const __m256i reg_nd = _mm256_set1_epi32(mint_type::Nd);
  do_butterfly<MOD, /*Inverse = */ false>(c1, n);
  do_butterfly<MOD, /*Inverse = */ false>(c2, n);
  for (int i = 0; i < n; i += 16) {
    const __m256i v0 = _mm256_load_si256((__m256i*)&c1[i]);
    const __m256i v1 = _mm256_load_si256((__m256i*)&c1[i + 8]);
    const __m256i u0 = _mm256_load_si256((__m256i*)&c2[i]);
    const __m256i u1 = _mm256_load_si256((__m256i*)&c2[i + 8]);
    const __m256i prod0 = mulmod256x32(v0, u0, reg_md, reg_nd);
    const __m256i prod1 = mulmod256x32(v1, u1, reg_md, reg_nd);
    _mm256_store_si256((__m256i*)&c1[i], prod0);
    _mm256_store_si256((__m256i*)&c1[i + 8], prod1);
    _mm256_store_si256((__m256i*)&c2[i], _mm256_setzero_si256());
    _mm256_store_si256((__m256i*)&c2[i + 8], _mm256_setzero_si256());
  }
  do_butterfly<MOD, /*Inverse = */ true>(c1, n);
}
// do garner's method
// precondition: the values of c1,c2 and c3 are already reduced
// precondition: n is power of 2 and sufficiently large
template <uint32_t MOD>
void consolidate(const uint32_t* c1, const uint32_t* c2, const uint32_t* c3,
                 uint32_t* out, int n) {
  constexpr uint32_t M1 = 167772161;
  constexpr uint32_t M2 = 469762049;
  constexpr uint32_t M3 = 754974721;
  constexpr auto I12 = montgomery_modint<M2>(M1).inv();
  constexpr auto I13 = montgomery_modint<M3>(M1).inv();
  constexpr auto I23 = montgomery_modint<M3>(M2).inv();
  constexpr auto C1 = montgomery_modint<MOD>(1);
  constexpr auto C2 = C1 * M1;
  constexpr auto C3 = C2 * M2;
  const __m256i reg_i12 = _mm256_set1_epi32(I12.get());
  const __m256i reg_i13 = _mm256_set1_epi32(I13.get());
  const __m256i reg_i23 = _mm256_set1_epi32(I23.get());
  const __m256i reg_co1 = _mm256_set1_epi32(C1.get());
  const __m256i reg_co2 = _mm256_set1_epi32(C2.get());
  const __m256i reg_co3 = _mm256_set1_epi32(C3.get());
  const __m256i reg_m2 = _mm256_set1_epi32(M2);
  const __m256i reg_m3 = _mm256_set1_epi32(M3);
  const __m256i reg_md = _mm256_set1_epi32(MOD);
  const __m256i reg_m2inv = _mm256_set1_epi32(montgomery_modint<M2>::Nd);
  const __m256i reg_m3inv = _mm256_set1_epi32(montgomery_modint<M3>::Nd);
  const __m256i reg_nd = _mm256_set1_epi32(montgomery_modint<MOD>::Nd);
  for (int i = 0; i < n; i += 16) {
    __m256i x1 = _mm256_load_si256((__m256i*)&c1[i]);
    __m256i y1 = _mm256_load_si256((__m256i*)&c1[i + 8]);
    __m256i x2 = _mm256_load_si256((__m256i*)&c2[i]);
    __m256i y2 = _mm256_load_si256((__m256i*)&c2[i + 8]);
    __m256i x3 = _mm256_load_si256((__m256i*)&c3[i]);
    __m256i y3 = _mm256_load_si256((__m256i*)&c3[i + 8]);
    x2 = mulmod256x32(submod256x32(x2, x1, reg_m2), reg_i12, reg_m2, reg_m2inv);
    y2 = mulmod256x32(submod256x32(y2, y1, reg_m2), reg_i12, reg_m2, reg_m2inv);
    const __m256i tmp1 =
        mulmod256x32(submod256x32(x3, x1, reg_m3), reg_i13, reg_m3, reg_m3inv);
    const __m256i tmp2 =
        mulmod256x32(submod256x32(y3, y1, reg_m3), reg_i13, reg_m3, reg_m3inv);
    x3 = mulmod256x32(submod256x32(tmp1, x2, reg_m3), reg_i23, reg_m3,
                      reg_m3inv);
    y3 = mulmod256x32(submod256x32(tmp2, y2, reg_m3), reg_i23, reg_m3,
                      reg_m3inv);
    x1 = mulmod256x32(x1, reg_co1, reg_md, reg_nd);
    y1 = mulmod256x32(y1, reg_co1, reg_md, reg_nd);
    x2 = mulmod256x32(x2, reg_co2, reg_md, reg_nd);
    y2 = mulmod256x32(y2, reg_co2, reg_md, reg_nd);
    x3 = mulmod256x32(x3, reg_co3, reg_md, reg_nd);
    y3 = mulmod256x32(y3, reg_co3, reg_md, reg_nd);
    x2 = addmod256x32(x1, x2, reg_md);
    y2 = addmod256x32(y1, y2, reg_md);
    x3 = addmod256x32(x2, x3, reg_md);
    y3 = addmod256x32(y2, y3, reg_md);
    _mm256_storeu_si256((__m256i*)&out[i], x3);
    _mm256_storeu_si256((__m256i*)&out[i + 8], y3);
  }
}
}  // namespace algo::detail::fft