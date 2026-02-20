#pragma once
#include <cassert>
#include <utility>
#include <vector>
// modules
#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"
/*
@internal/base/bit-base.hpp
@internal/base/typing.hpp
@internal/base/def.hpp
*/
// makecode
namespace algo {
/**
 * \brief Fenwick tree implementation. It is intended for prefix query
 * (such as prefix sum,prefix xor etc).
 *
 * \tparam Monoid Value type
 * \tparam BinOp Binary operation
 * \tparam Id Identity element
 */
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
class fenwick_tree {
public:
  fenwick_tree() noexcept = default;
  fenwick_tree(int n, BinOp binop, Id id) noexcept
      : tree_(n + 1, id()), binop_{std::move(binop)}, id_{std::move(id)} {}
  fenwick_tree(int n, const Monoid& def, BinOp binop, Id id) noexcept
      : tree_(n + 1, def), binop_{std::move(binop)}, id_{std::move(id)} {
    tree_[0] = id_();
    for (int i = 1; i <= n; ++i) {
      const int j = i + detail::blsi(i);
      if (j <= n) { tree_[j] = binop_(tree_[i], tree_[j]); }
    }
  }
  /**
   * \brief Update value at specified position.Let x be the value before
   * the update, and y be the new value, then x op z = y should hold,
   * otherwise the update is not valid anymore. (e.g. for addition query,
   * z+x=y => z=y-x)
   *
   * \param pos Position where the update is performed
   * \param z z op x = y (where op is BinOp) should hold
   */
  void apply(int pos, const Monoid& z) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    for (int i = pos + 1; i <= size(); i += detail::blsi(i)) {
      tree_[i] = binop_(tree_[i], z);
    }
  }
  /**
   * \brief Return the product F(a[0],a[1],...a[r-1])
   *
   * \param r upper bound (exclusive)
   * \return Product of F
   */
  Monoid query(int r) const {
#if !defined(NDEBUG)
    assert(0 <= r && r <= size());
#endif
    Monoid res = id_();
    for (int i = r; i > 0; i = detail::blsr(i)) { res = binop_(tree_[i], res); }
    return res;
  }
  int size() const { return static_cast<int>(tree_.size()) - 1; }
private:
  // The whole structure is 1 indexed
  std::vector<Monoid> tree_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
};
}  // namespace algo