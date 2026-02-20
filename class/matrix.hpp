#pragma once

#include <array>
#include <cassert>
#include <concepts>
#include <cstdlib>

#include "../internal/base/numeric-base.hpp"
#include "../internal/base/typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/numeric-base.hpp
*/
// makecode
namespace algo {
template <std::semiregular Ring, int Nd = 2>
struct basic_matrix : std::array<std::array<Ring, Nd>, Nd> {
private:
  static_assert(0 <= Nd, "Invalid dimention");
  using proxy_type = std::array<Ring, Nd>;
  using base_type = std::array<proxy_type, Nd>;
public:
  constexpr basic_matrix& operator+=(const basic_matrix& rhs) {
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { (*this)[i][j] += rhs[i][j]; }
    }
    return *this;
  }
  constexpr basic_matrix& operator-=(const basic_matrix& rhs) {
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { (*this)[i][j] -= rhs[i][j]; }
    }
    return *this;
  }
  constexpr basic_matrix& operator*=(const basic_matrix& rhs) {
    for (int i = 0; i < Nd; ++i) {
      proxy_type dest{};
      for (int k = 0; k < Nd; ++k) {
        for (int j = 0; j < Nd; ++j) { dest[j] += (*this)[i][k] * rhs[k][j]; }
      }
      (*this)[i] = dest;
    }
    return *this;
  }
  constexpr basic_matrix& operator*=(Ring rhs) {
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { (*this)[i][j] *= rhs; }
    }
    return *this;
  }
  constexpr basic_matrix operator-() const {
    basic_matrix res{};
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { res[i][j] = -(*this)[i][j]; }
    }
    return res;
  }
  constexpr basic_matrix transpose() const {
    basic_matrix res{};
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { res[i][j] = (*this)[j][i]; }
    }
    return res;
  }
  constexpr Ring trace() const {
    Ring res{};
    for (int i = 0; i < Nd; ++i) { res += (*this)[i][i]; }
    return res;
  }
  template <detail::integer Tp>
  constexpr basic_matrix pow(Tp n) const {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    return detail::binary_exp(*this, n, identity());
  }
  constexpr friend auto operator+(basic_matrix lhs, const basic_matrix& rhs)
      -> basic_matrix {
    lhs += rhs;
    return lhs;
  }
  constexpr friend auto operator-(basic_matrix lhs, const basic_matrix& rhs)
      -> basic_matrix {
    lhs -= rhs;
    return lhs;
  }
  constexpr friend auto operator*(basic_matrix lhs, const basic_matrix& rhs)
      -> basic_matrix {
    lhs *= rhs;
    return lhs;
  }
  constexpr friend auto operator*(basic_matrix lhs, Ring rhs) -> basic_matrix {
    lhs *= rhs;
    return lhs;
  }
  constexpr friend auto operator*(Ring lhs, basic_matrix rhs) -> basic_matrix {
    rhs *= lhs;
    return rhs;
  }
  static constexpr auto identity() -> basic_matrix {
    basic_matrix res{};
    for (int i = 0; i < Nd; ++i) { res[i][i] = Ring(1); }
    return res;
  }
  static constexpr auto zeros() -> basic_matrix { return {}; }
  static constexpr auto ones() -> basic_matrix {
    basic_matrix res{};
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { res[i][j] = Ring(1); }
    }
    return res;
  }
};
template <std::semiregular Ring, int Nd = 2>
struct basic_vector : std::array<Ring, Nd> {
private:
  using base_type = std::array<Ring, Nd>;
public:
  constexpr basic_vector& operator+=(const basic_vector& rhs) {
    for (int i = 0; i < Nd; ++i) { (*this)[i] += rhs[i]; }
    return *this;
  }
  constexpr basic_vector& operator-=(const basic_vector& rhs) {
    for (int i = 0; i < Nd; ++i) { (*this)[i] -= rhs[i]; }
    return *this;
  }
  constexpr basic_vector& operator*=(Ring rhs) {
    for (int i = 0; i < Nd; ++i) { (*this)[i] *= rhs; }
    return *this;
  }
  constexpr basic_vector operator-() const {
    basic_vector res{};
    for (int i = 0; i < Nd; ++i) { res[i] = -(*this)[i]; }
    return res;
  }
  constexpr friend auto operator+(basic_vector lhs, const basic_vector& rhs)
      -> basic_vector {
    lhs += rhs;
    return lhs;
  }
  constexpr friend auto operator-(basic_vector lhs, const basic_vector& rhs)
      -> basic_vector {
    lhs -= rhs;
    return lhs;
  }
  constexpr friend auto operator*(basic_vector lhs, Ring rhs) -> basic_vector {
    lhs *= rhs;
    return lhs;
  }
  constexpr friend auto operator*(Ring lhs, basic_vector rhs) -> basic_vector {
    rhs *= lhs;
    return rhs;
  }
  constexpr friend auto operator*(const basic_matrix<Ring, Nd>& lhs,
                                  const basic_vector& rhs) -> basic_vector {
    basic_vector res{};
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) { res[i] += lhs[i][j] * rhs[j]; }
    }
    return res;
  }
};
}  // namespace algo