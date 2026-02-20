#pragma once
#include <cassert>
#include <utility>
#include <vector>

#include "../internal/base/def.hpp"
/*
@internal/base/def.hpp
*/

// makecode
namespace algo {
/**
 * \brief Disjoint Set
 *
 */
class dsu {
public:
  dsu() noexcept = default;
  explicit dsu(int n) noexcept : parent_(n, -1) {}
  bool merge(int a, int b) {
#if !defined(NDEBUG)
    assert(0 <= a && a < num_nodes());
    assert(0 <= b && b < num_nodes());
#endif
    int sa = find(a), sb = find(b);
    if (sa != sb) {
      if (-parent_[sa] < -parent_[sb]) { std::swap(sa, sb); }
      parent_[sa] += parent_[sb];
      parent_[sb] = sa;
      return true;
    }
    return false;
  }
  int find(int a) {
#if !defined(NDEBUG)
    assert(0 <= a && a < num_nodes());
#endif
    if (parent_[a] >= 0) {
      parent_[a] = find(parent_[a]);
      return parent_[a];
    }
    return a;
  }
  int size(int a) {
#if !defined(NDEBUG)
    assert(0 <= a && a < num_nodes());
#endif
    return -parent_[find(a)];
  }
  bool same(int a, int b) {
#if !defined(NDEBUG)
    assert(0 <= a && a < num_nodes());
    assert(0 <= b && b < num_nodes());
#endif
    return find(a) == find(b);
  }
  int num_nodes() const { return static_cast<int>(parent_.size()); }
private:
  std::vector<int> parent_;
};
}  // namespace algo