#pragma once
#include <algorithm>
#include <cassert>
#include <iterator>
#include <vector>
// modules
#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
#include "./grid2d.hpp"

/*
@internal/base/bit-base.hpp
@internal/base/typing.hpp
@internal/base/def.hpp
@class/grid2d.hpp
*/
// makecode
namespace algo {
/**
 * \brief typical implementation of sparse table.
Here i assert followings hold:
1.F(x, y, z) = F(x, F(y, z)) = F(F(x, y), z)
2.F(x, y) = F(F(x, x), y) = F(x, F(y, y)) = F(y, x)
 *
 * \tparam SemiGroup value type
 * \tparam BinOp Binary operation function
 */
template <std::semiregular SemiGroup,
          detail::function<SemiGroup(SemiGroup, SemiGroup)> BinOp>
class sparse_table {
public:
  sparse_table() noexcept = default;
  sparse_table(int n, const SemiGroup &def, BinOp binop) noexcept
      : table_(detail::bit_width(n), n, def), binop_{std::move(binop)} {
    std::fill_n(table_.data(), n, def);
    calc_table();
  }
  template <detail::input_iterator<SemiGroup> IIter>
  sparse_table(IIter first, IIter last, BinOp binop) noexcept
      : table_(detail::bit_width(std::distance(first, last)),
               static_cast<int>(std::distance(first, last))),
        binop_{std::move(binop)} {
    std::copy(first, last, table_.data());
    calc_table();
  }
  template <std::input_iterator IIter,
            detail::function<SemiGroup(std::iter_value_t<IIter>)> Gen>
  sparse_table(IIter first, IIter last, BinOp binop, Gen gen) noexcept
      : table_(detail::bit_width(std::distance(first, last)),
               static_cast<int>(std::distance(first, last))),
        binop_{std::move(binop)} {
    std::transform(first, last, table_.data(), std::move(gen));
    calc_table();
  }
  /**
   * \brief Return the product of BinOp over elements in range [l,r) in
   * O(1) assuming the binary operation is idempotent.
   *
   * \param l Left boundary
   * \param r Right boundary (exclusive)
   * \return SemiGroup
   */
  SemiGroup query(int l, int r) const {
#if !defined(NDEBUG)
    assert(0 <= l && l < r && r <= size());
#endif
    const int lg = detail::floor_log2(r - l);
    return binop_(table_[lg][l], table_[lg][r - (1 << lg)]);
  }
  /**
   * \brief Return the product of BinOp over elements in range [l,r) in
   * O(log(r - l)).The binary operation may not be idempotent
   *
   * \param l
   * \param r
   * \return
   */
  SemiGroup query_incremental(int l, int r) const {
#if !defined(NDEBUG)
    assert(l < r && l >= 0 && r <= size());
#endif
    const int istep = detail::floor_log2(r - l);
    SemiGroup res = table_[istep][l];
    int base = l + (1 << istep), len = r - base;
    while (len > 0) {
      const int nstep = detail::floor_log2(len);
      res = binop_(res, table_[nstep][base]);
      base += 1 << nstep;
      len -= 1 << nstep;
    }
    return res;
  }
  int size() const { return table_.col_size(); }
private:
  void calc_table() {
    for (int j = 1; j < table_.row_size(); ++j) {
      for (int i = 0; i + (1 << j) <= table_.col_size(); ++i) {
        table_[j][i] =
            binop_(table_[j - 1][i], table_[j - 1][i + (1 << (j - 1))]);
      }
    }
  }
  algo::grid2d<SemiGroup> table_;
  [[no_unique_address]] BinOp binop_;
};
}  // namespace algo