#pragma once

#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <utility>
#include <vector>
// modules
#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/bit-base.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo {
/**
 * \brief Typical segment tree implementation for range query and point
 * update.
 *
 * \tparam Monoid
 * \tparam BinOp
 * \tparam Id
 */
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
class segtree {
public:
  segtree() noexcept = default;
  segtree(int n, BinOp binop, Id id) noexcept
      : n_{n},
        lg2_{detail::ceil_log2(n)},
        size_{1 << lg2_},
        tree_(2 * size_, id()),
        binop_{std::move(binop)},
        id_{std::move(id)} {}
  segtree(int n, const Monoid& init, BinOp binop, Id id) noexcept
      : segtree(n, std::move(binop), std::move(id)) {
    std::fill_n(tree_.begin() + size_, size(), init);
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  template <detail::input_iterator<Monoid> IIter>
  segtree(IIter first, IIter last, BinOp binop, Id id) noexcept
      : segtree(static_cast<int>(std::distance(first, last)), std::move(binop),
                std::move(id)) {
    std::copy(first, last, tree_.begin() + size_);
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  template <std::input_iterator IIter,
            detail::function<Monoid(std::iter_value_t<IIter>)> Gen>
  segtree(IIter first, IIter last, BinOp binop, Id id, Gen gen) noexcept
      : segtree(static_cast<int>(std::distance(first, last)), std::move(binop),
                std::move(id)) {
    std::transform(first, last, tree_.begin() + size_, std::move(gen));
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  /**
   * \brief Set t[pos] to new_val
   *
   * \param pos 0<=pos<n_ should be satisfied
   * \param new_val The value to set
   */
  void set(int pos, Monoid new_val) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    tree_[pos] = std::move(new_val);
    for (int i = 1; i <= lg2_; ++i) { update(pos >> i); }
  }

  /**
   * \brief Return t[pos]
   *
   * \param pos
   * \return  value at pos
   */
  Monoid get(int pos) const {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    return tree_[pos + size_];
  }
  /**
   * \brief Set t[pos] to F(t[pos], z).
   *
   * \param pos
   * \param z
   */
  void apply(int pos, Monoid z) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    tree_[pos] = binop_(tree_[pos], std::move(z));
    for (int i = 1; i <= lg2_; ++i) { update(pos >> i); }
  }
  /**
   * \brief Return the product of BinOp(a[l],a[l+1],...,a[r-1])
   *
   * \param l Left side boundary
   * \param r Right side boundary (exclusive)
   * \return Product of BinOp
   */
  Monoid query(int l, int r) const {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    Monoid resl = id_(), resr = id_();
    l += size_ - 1;
    r += size_ + 0;
    const int mask = detail::floor_pow2(l ^ r) - 1;
    for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
      resl = binop_(resl, tree_[(l >> detail::count_tz(v)) + 1]);
    }
    for (int v = mask & r; v != 0; v = detail::blsr(v)) {
      resr = binop_(tree_[(r >> detail::count_tz(v)) - 1], resr);
    }
    return binop_(std::move(resl), std::move(resr));
  }
  /**
   * \brief Return all product of BinOp(a[0],a[1],...a[n-1])
   *
   * \return Product of BinOp
   */
  Monoid query() const { return tree_[1]; }
  /**
   * \brief Return the largest index r satisfying both of the followings:
   * 1.r=l or pred(f(a[l],a[l+1],...a[r-1]))=true
   * 2.r=n or pred(f[a[l],a[l+1],...a[r]])=false
   *
   * \tparam Pred some monotonic function returning bool value
   * \param l Left boundary
   * \param pred pred(e())=true should be hold
   * \return
   */
  template <std::predicate<Monoid> Pred>
  int bisect_right(int l, const Pred& pred) const {
#if !defined(NDEBUG)
    assert(0 <= l && l <= size());
    assert(pred(id_()));
#endif
    if (l < size()) [[likely]] {
      l += size_;
      Monoid cur = id_();
      do {
        l >>= detail::count_tz(l);
        if (!pred(binop_(cur, tree_[l]))) {
          while (l < size_) {
            l = 2 * l;
            if (Monoid tmp = binop_(cur, tree_[l]); pred(tmp)) {
              cur = std::move(tmp), ++l;
            }
          }
          return l - size_;
        }
        cur = binop_(cur, tree_[l++]);
      } while (detail::popcount(l) != 1);
    }
    return size();
  }
  /**
   * \brief Return smallest index l satisfying both of the followings:
   * 1.l=r or pred(f(a[l],a[l+1],...a[r-1]))=true
   * 2.l=0 or pred(f(a[l-1],a[l],...a[r-1]))=false
   *
   * \tparam Pred Unary predicate accepting Monoid
   * \param r Right side boundary
   * \param pred pred(e())=true should be hold
   * \return
   */
  template <std::predicate<Monoid> Pred>
  int bisect_left(int r, const Pred& pred) const {
#if !defined(NDEBUG)
    assert(0 <= r && r <= size());
#endif
    if (r > 0) [[likely]] {
      r += size_ - 1;
      Monoid cur = id_();
      do {
        r >>= detail::count_tz(~r);
        if (!pred(binop_(tree_[r], cur))) {
          while (r < size_) {
            r = 2 * r + 1;
            if (Monoid tmp = binop_(tree_[r], cur); pred(tmp)) {
              cur = std::move(tmp), --r;
            }
          }
          return r + 1 - size_;
        }
        cur = binop_(tree_[r], cur);
      } while (detail::popcount(r--) != 1);
    }
    return 0;
  }
  int size() const { return n_; }
private:
  void update(int v) { tree_[v] = binop_(tree_[2 * v], tree_[2 * v + 1]); }
  int n_;
  int lg2_;
  int size_;
  std::vector<Monoid> tree_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
};
}  // namespace algo