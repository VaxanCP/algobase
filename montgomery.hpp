#pragma once

#include "./class/montgomery.hpp"

/*
@class/montgomery.hpp
*/

#if defined(_DEBUG)
#  include <fmt/format.h>
template <uint32_t MOD>
struct fmt::formatter<algo::montgomery_modint<MOD>> : fmt::formatter<uint32_t> {
  constexpr auto format(algo::montgomery_modint<MOD> m,
                        fmt::format_context& ctx) const {
    return fmt::formatter<uint32_t>::format(m.value(), ctx);
  }
};
#endif