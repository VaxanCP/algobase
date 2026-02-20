#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <utility>
#include <vector>

#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
/*
@internal/base/typing.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo {
template <std::semiregular Tp>
class grid2d {
public:
  grid2d() noexcept = default;
  grid2d(int row, int col, const Tp& def = Tp{}) noexcept
      : row_{row}, col_{col}, data_(total_size(), def) {}
  std::span<const Tp> operator[](int i) const {
#if !defined(NDEBUG)
    assert(0 <= i && i < row_size());
#endif
    return {data_.cbegin() + i * col_size(),
            data_.cbegin() + (i + 1) * col_size()};
  }
  std::span<Tp> operator[](int i) {
#if !defined(NDEBUG)
    assert(0 <= i && i < row_size());
#endif
    return {data_.begin() + i * col_size(),
            data_.begin() + (i + 1) * col_size()};
  }
  const Tp& operator()(int i, int j) const {
#if !defined(NDEBUG)
    assert(0 <= i && i < row_size());
    assert(0 <= j && j < col_size());
#endif
    return data_[i * col_size() + j];
  }
  Tp& operator()(int i, int j) {
#if !defined(NDEBUG)
    assert(0 <= i && i < row_size());
    assert(0 <= j && j < col_size());
#endif
    return data_[i * col_size() + j];
  }
  int total_size() const { return row_ * col_; }
  int row_size() const { return row_; }
  int col_size() const { return col_; }
  const Tp* data() const { return data_.data(); }
  Tp* data() { return data_.data(); }
private:
  int row_;
  int col_;
  std::vector<Tp> data_;
};
}  // namespace algo

#if defined(_DEBUG)
#  include <fmt/format.h>
#  include <fmt/ranges.h>
template <typename Tp>
struct fmt::formatter<algo::grid2d<Tp>> {
  auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }
  template <typename FormatContext>
  auto format(const algo::grid2d<Tp>& v, FormatContext& ctx) const {
    auto buf = fmt::memory_buffer();
    fmt::format_to(std::back_inserter(buf), "[\n");
    for (int i = 0; i < v.row_size(); ++i) {
      fmt::format_to(std::back_inserter(buf), "{}", v[i]);
      fmt::format_to(std::back_inserter(buf), "\n");
    }
    fmt::format_to(std::back_inserter(buf), "]");
    return fmt::format_to(ctx.out(), "{}", fmt::to_string(buf));
  }
};
#endif
