#pragma once

#include <cassert>
#include <cmath>

#include "../internal/base/def.hpp"
#include "../internal/base/math-base.hpp"
#include "../internal/base/typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/math-base.hpp
@interanl/base/def.hpp
*/

// makecode
namespace algo {
template <detail::arithmetic Tp>
struct point2d {
  Tp x;
  Tp y;
  constexpr point2d& operator+=(point2d p) {
    x += p.x;
    y += p.y;
    return *this;
  }
  constexpr point2d& operator-=(point2d p) {
    x -= p.x;
    y -= p.y;
    return *this;
  }
  constexpr point2d& operator*=(Tp s) {
    x *= s;
    y *= s;
    return *this;
  }
  constexpr point2d& operator*=(point2d p) {
    x *= p.x;
    y *= p.y;
    return *this;
  }
  constexpr point2d& operator/=(Tp s) {
#if !defined(NDEBUG)
    assert(s != Tp(0));
#endif
    x /= s;
    y /= s;
    return *this;
  }
  constexpr point2d& operator/=(point2d p) {
#if !defined(NDEBUG)
    assert(p.x != Tp(0));
    assert(p.y != Tp(0));
#endif
    x /= p.x;
    y /= p.y;
    return *this;
  }
  friend constexpr auto operator-(point2d p) -> point2d {
    p.x = -p.x;
    p.y = -p.y;
    return p;
  }
  friend constexpr auto operator+(point2d lhs, point2d rhs) -> point2d {
    lhs += rhs;
    return lhs;
  }
  friend constexpr auto operator-(point2d lhs, point2d rhs) -> point2d {
    lhs -= rhs;
    return lhs;
  }
  friend constexpr auto operator*(point2d lhs, point2d rhs) -> point2d {
    lhs *= rhs;
    return lhs;
  }
  friend constexpr auto operator*(point2d lhs, Tp rhs) -> point2d {
    lhs *= rhs;
    return lhs;
  }
  friend constexpr auto operator*(Tp lhs, point2d rhs) -> point2d {
    rhs *= lhs;
    return rhs;
  }
  friend constexpr auto operator/(point2d lhs, point2d rhs) -> point2d {
    lhs /= rhs;
    return lhs;
  }
  friend constexpr auto operator/(point2d lhs, Tp rhs) -> point2d {
    lhs /= rhs;
    return lhs;
  }
  // Comparison operator
  auto operator<=>(const point2d&) const = default;
  // Dot product
  constexpr Tp dot(point2d p) const { return x * p.x + y * p.y; }
  constexpr Tp dot(point2d p1, point2d p2) const {
    const point2d p01 = p1 - *this;
    const point2d p02 = p2 - *this;
    return p01.dot(p02);
  }
  constexpr auto safe_dot(point2d p) const
    requires detail::integer<Tp>
  {
    using promoted_type = detail::imul_result_t<Tp>;
    const auto xx = promoted_type(x) * p.x;
    const auto yy = promoted_type(y) * p.y;
    return xx + yy;
  }
  constexpr auto safe_dot(point2d p1, point2d p2) const
    requires detail::integer<Tp>
  {
    const point2d p01 = p1 - *this;
    const point2d p02 = p2 - *this;
    return p01.safe_dot(p02);
  }
  constexpr Tp cross(point2d p) const { return x * p.y - y * p.x; }
  constexpr Tp cross(point2d p1, point2d p2) const {
    const point2d p01 = p1 - *this;
    const point2d p02 = p2 - *this;
    return p01.cross(p02);
  }
  constexpr auto safe_cross(point2d p) const
    requires detail::integer<Tp>
  {
    using promoted_type = detail::imul_result_t<Tp>;
    const auto x1y2 = promoted_type(x) * p.y;
    const auto y1x2 = promoted_type(y) * p.x;
    return x1y2 - y1x2;
  }
  constexpr auto safe_cross(point2d p1, point2d p2) const
    requires detail::integer<Tp>
  {
    const point2d p01 = p1 - *this;
    const point2d p02 = p2 - *this;
    return p01.safe_cross(p02);
  }
  constexpr Tp norm(point2d p) const { return dot(*this); }
  constexpr auto safe_norm() const
    requires detail::integer<Tp>
  {
    return safe_dot(*this);
  }
  constexpr point2d rotate_left() const { return {-y, x}; }
  constexpr point2d rotate_right() const { return {y, -x}; }
};
}  // namespace algo