#pragma once

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "../internal/base/bit-base.hpp"
#include "../internal/base/def.hpp"
#include "../internal/base/typing.hpp"

/*
@internal/base/def.hpp
@internal/base/bit-base.hpp
@internal/base/typing.hpp
*/

// makecode
namespace algo {
// Forward declaration
template <typename Monoid, typename Mact, typename Map, typename Comp,
          typename CompId, typename BinOp = detail::monostate,
          typename BinOpId = detail::monostate>
class lazy_segtree;

/**
 * \brief Segment tree implementation for both range update and range
 * query. To perform operations successfully, followings must hold: For
 * Monoid S equipped with a binary operation F, Monoid A with a binary
 * operation G, and a mapping function H : S x G -> S, then
 * 1. H(x,e) = x for x in S where e is the identity element of A.
 * 2. H(H(x,a),b) = H(x,G(a,b)) for x,y in S and a,b in A.
 * 3. F(H(x,a),H(y,a)) = H(F(x,y),a) for x,y in S and a in A.
 *
 * \tparam Monoid
 * \tparam Mact
 * \tparam Map
 * \tparam Comp
 * \tparam CompId
 * \tparam BinOp
 * \tparam BinOpId
 */
template <std::semiregular Monoid, std::semiregular Mact,
          detail::function<Monoid(Monoid, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> BinOpId>
class lazy_segtree<Monoid, Mact, Map, Comp, CompId, BinOp, BinOpId> {
public:
  lazy_segtree() noexcept = default;
  /**
   * \brief Construct tree of size n with passed functors
   *
   * \param n The size of tree
   * \param map The mapping function
   * \param comp The binary operation of A
   * \param comp_id The identity of A
   * \param binop The binary operation of S
   * \param binop_id The identity of S
   */
  lazy_segtree(int n, Map map, Comp comp, CompId comp_id, BinOp binop,
               BinOpId binop_id) noexcept
      : n_{n},
        lg2_{detail::ceil_log2(n)},
        size_{1 << lg2_},
        tree_(2 * size_, binop_id()),
        lazy_(size_, comp_id()),
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)},
        binop_{std::move(binop)},
        binop_id_{std::move(binop_id)} {}
  lazy_segtree(int n, const Monoid& def, Map map, Comp comp, CompId comp_id,
               BinOp binop, BinOpId binop_id) noexcept
      : lazy_segtree(n, std::move(map), std::move(comp), std::move(comp_id),
                     std::move(binop), std::move(binop_id)) {
    std::fill_n(tree_.begin() + size_, size(), def);
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  template <detail::input_iterator<Monoid> IIter>
  lazy_segtree(IIter first, IIter last, Map map, Comp comp, CompId comp_id,
               BinOp binop, BinOpId binop_id) noexcept
      : lazy_segtree(static_cast<int>(std::distance(first, last)),
                     std::move(map), std::move(comp), std::move(comp_id),
                     std::move(binop), std::move(binop_id)) {
    std::copy(first, last, tree_.begin() + size_);
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  template <std::input_iterator IIter,
            detail::function<Monoid(std::iter_value_t<IIter>)> Gen>
  lazy_segtree(IIter first, IIter last, Map map, Comp comp, CompId comp_id,
               BinOp binop, BinOpId binop_id, Gen gen) noexcept
      : lazy_segtree(static_cast<int>(std::distance(first, last)),
                     std::move(map), std::move(comp), std::move(comp_id),
                     std::move(binop), std::move(binop_id)) {
    std::transform(first, last, tree_.begin() + size_, std::move(gen));
    for (int i = size_ - 1; i > 0; --i) { update(i); }
  }
  /**
   * \brief Set \p new_val to element at position \p pos
   *
   * \param pos Position of the element.
   * \param new_val New value to set.
   */
  void set(int pos, Monoid new_val) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = lg2_; i > 0; --i) { push(pos >> i); }
    tree_[pos] = std::move(new_val);
    for (int i = 1; i <= lg2_; ++i) { update(pos >> i); }
  }
  /**
   * \brief Get element at specified position.
   *
   * \param pos Position of the element.
   * \return Value at pos
   */
  Monoid get(int pos) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = lg2_; i > 0; --i) { push(pos >> i); }
    return tree_[pos];
  }
  /**
   * \brief Apply non-commutative action to element at specified position
   *
   * \param pos
   * \param arg
   */
  void apply(int pos, Mact arg) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = lg2_; i > 0; --i) { push(pos >> i); }
    tree_[pos] = map_(tree_[pos], std::move(arg));
    for (int i = 1; i <= lg2_; ++i) { update(pos >> i); }
  }
  /**
   * \brief Apply commutative action to element at specified position
   *
   * \param pos
   * \param arg
   */
  void apply_commutative(int pos, Mact arg) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    tree_[pos] = map_(tree_[pos], std::move(arg));
    for (int i = 1; i <= lg2_; ++i) { push_update(pos >> i); }
  }
  /**
   * \brief Apply non-commutative action for elements in range [l, r)
   *
   * \param l
   * \param r
   * \param arg
   */
  void apply(int l, int r, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) {
      l += size_ - 1;
      r += size_ + 0;
      const int tzl = detail::count_tz(~l);
      const int tzr = detail::count_tz(r);
      const int mask = detail::floor_pow2(l ^ r) - 1;
      for (int i = lg2_; i > tzl; --i) { push(l >> i); }
      for (int i = lg2_; i > tzr; --i) { push(r >> i); }
      for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
        propagate((l >> detail::count_tz(v)) + 1, arg);
      }
      for (int v = mask & r; v != 0; v = detail::blsr(v)) {
        propagate((r >> detail::count_tz(v)) - 1, arg);
      }
      for (int i = tzl + 1; i <= lg2_; ++i) { update(l >> i); }
      for (int i = tzr + 1; i <= lg2_; ++i) { update(r >> i); }
    }
  }
  /**
   * \brief Apply commutative action for elements in range [l,r).
   *
   * \param l A left bound.
   * \param r A right bound (exclusive)
   * \param arg An argument which applied to the elements
   */
  void apply_commutative(int l, int r, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) {
      l += size_ - 1;
      r += size_ + 0;
      const int tzl = detail::count_tz(~l);
      const int tzr = detail::count_tz(r);
      const int mask = detail::floor_pow2(l ^ r) - 1;
      for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
        propagate((l >> detail::count_tz(v)) + 1, arg);
      }
      for (int v = mask & r; v != 0; v = detail::blsr(v)) {
        propagate((r >> detail::count_tz(v)) - 1, arg);
      }
      for (int i = tzl + 1; i <= lg2_; ++i) { push_update(l >> i); }
      for (int i = tzr + 1; i <= lg2_; ++i) { push_update(r >> i); }
    }
  }
  /**
   * \brief Compute product of the elements in range [l, r). If l = r,
   * return identity.
   *
   * \param l A left bound
   * \param r A right bound (exclusive)
   * \return Monoid The product of the values of elements in range [l, r)
   */
  Monoid query(int l, int r) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) [[likely]] {
      Monoid resl = binop_id_(), resr = binop_id_();
      l += size_ - 1;
      r += size_ + 0;
      const int tzl = detail::count_tz(~l);
      const int tzr = detail::count_tz(r);
      const int mask = detail::floor_pow2(l ^ r) - 1;
      for (int i = lg2_; i > tzl; --i) { push(l >> i); }
      for (int i = lg2_; i > tzr; --i) { push(r >> i); }
      for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
        resl = binop_(resl, tree_[(l >> detail::count_tz(v)) + 1]);
      }
      for (int v = mask & r; v != 0; v = detail::blsr(v)) {
        resr = binop_(tree_[(r >> detail::count_tz(v)) - 1], resr);
      }
      return binop_(std::move(resl), std::move(resr));
    }
    return binop_id_();
  }
  /**
   * \brief Compute the product of all elements.
   *
   * \return The product of all elements.
   */
  Monoid query() const { return tree_[1]; }
  /**
   * \brief
   *
   * \tparam Pred
   * \param l
   * \param pred
   * \return
   */
  template <std::predicate<Monoid> Pred>
  int bisect_right(int l, const Pred& pred) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= size());
    assert(pred(binop_id_()));
#endif
    if (l < size()) [[likely]] {
      l += size_;
      for (int i = lg2_; i > 0; --i) { push(l >> i); }
      Monoid cur = binop_id_();
      do {
        l >>= detail::count_tz(l);
        if (!pred(binop_(cur, tree_[l]))) {
          while (l < size_) {
            push(l);
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
   * \brief
   *
   * \tparam Pred
   * \param r
   * \param pred
   * \return
   */
  template <std::predicate<Monoid> Pred>
  int bisect_left(int r, const Pred& pred) {
#if !defined(NDEBUG)
    assert(0 <= r && r <= size());
    assert(pred(binop_id_()));
#endif
    if (r > 0) [[likely]] {
      r += size_ - 1;
      for (int i = lg2_; i > 0; --i) { push(r >> i); }
      Monoid cur = binop_id_();
      do {
        r >>= detail::count_tz(~r);
        if (!pred(binop_(tree_[r], cur))) {
          while (r < size_) {
            push(r);
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
  void push_update(int v) {
    update(v);
    tree_[v] = map_(tree_[v], lazy_[v]);
  }
  void propagate(int v, const Mact& arg) {
    tree_[v] = map_(tree_[v], arg);
    if (v < size_) { lazy_[v] = comp_(lazy_[v], arg); }
  }
  void push(int v) {
    propagate(2 * v, lazy_[v]);
    propagate(2 * v + 1, lazy_[v]);
    lazy_[v] = comp_id_();
  }
  int n_;     // size of the input
  int lg2_;   // [log(n)]
  int size_;  // size of the underlying
  std::vector<Monoid> tree_;
  std::vector<Mact> lazy_;
  [[no_unique_address]] Map map_;
  [[no_unique_address]] Comp comp_;
  [[no_unique_address]] CompId comp_id_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] BinOpId binop_id_;
};

/**
 * \brief Segment tree implementation for range update and point query.
 *
 * \tparam Tp
 * \tparam Mact
 * \tparam Map
 * \tparam Comp
 * \tparam CompId
 */
template <std::semiregular Tp, std::semiregular Mact,
          detail::function<Tp(Tp, Mact)> Map,
          detail::function<Mact(Mact, Mact)> Comp,
          detail::function<Mact(void)> CompId>
class lazy_segtree<Tp, Mact, Map, Comp, CompId, detail::monostate,
                   detail::monostate> {
public:
  lazy_segtree() noexcept = default;
  lazy_segtree(int n, Map map, Comp comp, CompId comp_id) noexcept
      : size_{n},
        tree_(size_),
        lazy_(size_, comp_id()),
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)} {}
  lazy_segtree(int n, const Tp& def, Map map, Comp comp,
               CompId comp_id) noexcept
      : size_{n},
        tree_(size_, def),
        lazy_(size_, comp_id()),
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)} {}
  lazy_segtree(std::vector<Tp> vec, Map map, Comp comp, CompId comp_id) noexcept
      : size_{static_cast<int>(vec.size())},
        tree_(std::move(vec)),
        lazy_(size_, comp_id()),
        map_{std::move(map)},
        comp_{std::move(comp)},
        comp_id_{std::move(comp_id)} {}
  template <detail::input_iterator<Tp> IIter>
  lazy_segtree(IIter first, IIter last, Map map, Comp comp,
               CompId comp_id) noexcept
      : lazy_segtree(static_cast<int>(std::distance(first, last)),
                     std::move(map), std::move(comp), std::move(comp_id)) {
    std::copy(first, last, tree_.begin());
  }
  template <std::input_iterator IIter,
            detail::function<Tp(std::iter_value_t<IIter>)> Gen>
  lazy_segtree(IIter first, IIter last, Map map, Comp comp, CompId comp_id,
               Gen gen) noexcept
      : lazy_segtree(static_cast<int>(std::distance(first, last)),
                     std::move(map), std::move(comp), std::move(comp_id)) {
    std::transform(first, last, tree_.begin(), std::move(gen));
  }
  void set(int pos, Tp new_val) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = detail::floor_log2(pos); i > 0; --i) { push(pos >> i); }
    tree_[pos - size_] = std::move(new_val);
  }
  Tp get(int pos) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = detail::floor_log2(pos); i > 0; --i) { push(pos >> i); }
    return tree_[pos - size_];
  }
  void apply(int pos, Mact arg) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    pos += size_;
    for (int i = detail::floor_log2(pos); i > 0; --i) { push(pos >> i); }
    tree_[pos - size_] = map_(tree_[pos - size_], std::move(arg));
  }
  void apply_commutative(int pos, Mact arg) {
#if !defined(NDEBUG)
    assert(0 <= pos && pos < size());
#endif
    tree_[pos] = map_(tree_[pos], std::move(arg));
  }
  void apply(int l, int r, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) {
      l += size_ - 1;
      r += size_ + 0;
      const int tzl = detail::count_tz(~l);
      const int tzr = detail::count_tz(r);
      const int mask = detail::floor_pow2(l ^ r) - 1;
      for (int i = detail::floor_log2(l); i > tzl; --i) { push(l >> i); }
      for (int i = detail::floor_log2(r); i > tzr; --i) { push(r >> i); }
      for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
        propagate((l >> detail::count_tz(v)) + 1, arg);
      }
      for (int v = mask & r; v != 0; v = detail::blsr(v)) {
        propagate((r >> detail::count_tz(v)) - 1, arg);
      }
    }
  }
  void apply_commutative(int l, int r, const Mact& arg) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) {
      l += size_ - 1;
      r += size_ + 0;
      const int mask = detail::floor_pow2(l ^ r) - 1;
      for (int v = mask & ~l; v != 0; v = detail::blsr(v)) {
        propagate((l >> detail::count_tz(v)) + 1, arg);
      }
      for (int v = mask & r; v != 0; v = detail::blsr(v)) {
        propagate((r >> detail::count_tz(v)) - 1, arg);
      }
    }
  }
  int size() const { return size_; }
private:
  void propagate(int v, Mact arg) {
    if (v < size_) [[likely]] {
      lazy_[v] = comp_(lazy_[v], std::move(arg));
    } else {
      tree_[v - size_] = map_(tree_[v - size_], std::move(arg));
    }
  }
  void push(int v) {
    propagate(2 * v, lazy_[v]);
    propagate(2 * v + 1, lazy_[v]);
    lazy_[v] = comp_id_();
  }
  int size_;
  std::vector<Tp> tree_;
  std::vector<Mact> lazy_;
  [[no_unique_address]] Map map_;
  [[no_unique_address]] Comp comp_;
  [[no_unique_address]] CompId comp_id_;
};
}  // namespace algo