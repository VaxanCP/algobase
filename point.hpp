#pragma once

#include "./class/point.hpp"

/*
@class/point.hpp
*/

#if defined(_DEBUG)
#  include <fmt/format.h>
template <typename Tp>
struct fmt::formatter<algo::point2d<Tp>> {
  constexpr auto parse(fmt::format_parse_context& cxt) { return cxt.begin(); }
  template <typename FormatContext>
  constexpr auto format(algo::point2d<Tp> p, FormatContext& ctx) {
    return fmt::format_to(ctx.out(), "({0},{1})", p.x, p.y);
  }
};
#endif