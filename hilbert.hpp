#pragma once

#include "./internal/base/bit-base.hpp"
#include "./internal/base/typing.hpp"

/*
@internal/base/bit-base.hpp
@internal/base/typing.hpp
*/

// makecode
namespace algo {
class hilbert_plane {
public:
  constexpr hilbert_plane() noexcept = default;
  constexpr hilbert_plane(int n) noexcept
      : lg2_{detail::ceil_log2(n)}, dim_{1 << lg2_} {}
  constexpr int64_t order(int x, int y) const {
#if !defined(NDEBUG)
    assert(0 <= x && x < dimension());
    assert(0 <= y && y < dimension());
#endif
    int64_t res = 0;
    for (int i = log_dim(); i-- > 0;) {
      const bool rx = (x >> i) & 1;
      const bool ry = (y >> i) & 1;
      res += (static_cast<int64_t>(1) << (2 * i)) * Rot[rx][ry];
      if (!rx) {
        x = ry ? dimension() - 1 - x : x;
        y = ry ? dimension() - 1 - y : y;
        std::swap(x, y);
      }
    }
    return res;
  }
  constexpr int dimension() const { return dim_; }
  constexpr int log_dim() const { return lg2_; }
private:
  static constexpr int Rot[2][2]{{0, 3}, {1, 2}};
  int lg2_;
  int dim_;
};
}  // namespace algo