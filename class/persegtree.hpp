#pragma once

#include <cassert>
#include <memory_resource>
#include <stack>
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
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
class persistent_segtree {
private:
  struct pseg_node {
    Monoid value;
    pseg_node* left = nullptr;
    pseg_node* right = nullptr;
  };
  struct pseg_node_allocator {
    using node_type = pseg_node;
    using node_pointer = node_type*;
    using const_node_pointer = const node_type*;
    using pool_resource_type = std::pmr::monotonic_buffer_resource;
    static node_pointer allocate(Monoid value) {
      void* mem = do_allocate();
      return ::new (mem) node_type(std::move(value));
    }
    static node_pointer clone(const_node_pointer source) {
      void* mem = do_allocate();
#if !defined(NDEBUG)
      assert(source != nullptr);
#endif
      return ::new (mem) node_type(*source);
    }
    static void* do_allocate() {
      return pool_.allocate(sizeof(node_type), alignof(node_type));
    }
  private:
    static inline pool_resource_type pool_;
  };
  using node_allocator = pseg_node_allocator;
  using node_type = typename node_allocator::node_type;
  using node_pointer = typename node_allocator::node_pointer;
  using const_node_pointer = typename node_allocator::const_node_pointer;
public:
  persistent_segtree() noexcept = default;
  persistent_segtree(int n, BinOp binop, Id id) noexcept
      : size_{n}, binop_{std::move(binop)}, id_{std::move(id)}, tree_{nullptr} {
    if (size_ != 0) { tree_.front() = build(n, id_()); }
  }
  persistent_segtree(int n, const Monoid& def, BinOp binop, Id id) noexcept
      : size_{n}, binop_{std::move(binop)}, id_{std::move(id)}, tree_{nullptr} {
    if (size_ != 0) { tree_.front() = build(n, def); }
  }
  template <detail::input_iterator<Monoid> IIter>
  persistent_segtree(IIter first, IIter last, BinOp binop, Id id) noexcept
      : size_{static_cast<int>(std::distance(first, last))},
        binop_{std::move(binop)},
        id_{std::move(id)},
        tree_{nullptr} {
    if (size_ != 0) { tree_.front() = build(first, size_); }
  }
  Monoid get(int pos, int version) const {
#if !defined(NDEBUG)
    assert(0 <= version && version <= latest());
    assert(0 <= pos && pos < size());
#endif
    return search(tree_[version], pos)->value;
  }
  Monoid get(int pos) const { return get(pos, latest()); }
  void set(int pos, Monoid new_val, int version) {
#if !defined(NDEBUG)
    assert(0 <= version && version <= latest());
    assert(0 <= pos && pos < size());
#endif
    node_pointer new_root = search_update(
        tree_[version], pos,
        [new_val = std::move(new_val)](Monoid) { return std::move(new_val); });
    add_latest(new_root);
  }
  void set(int pos, Monoid new_val) { set(pos, std::move(new_val), latest()); }
  void apply(int pos, Monoid z, int version) {
#if !defined(NDEBUG)
    assert(0 <= version && version <= latest());
    assert(0 <= pos && pos < size());
#endif
    node_pointer new_root = search_update(
        tree_[version], pos, [z = std::move(z), this](Monoid val) {
          return binop_(std::move(val), std::move(z));
        });
    add_latest(new_root);
  }
  void apply(int pos, Monoid z) { apply(pos, std::move(z), latest()); }
  Monoid query(int l, int r, int version) const {
#if !defined(NDEBUG)
    assert(0 <= version && version <= latest());
    assert(0 <= l && l <= r && r <= size());
#endif
    if (l < r) [[likely]] {
      return product(tree_[version], 0, size(), l, r);
    } else {
      return id_();
    }
  }
  Monoid query(int l, int r) const { return query(l, r, latest()); }
  Monoid query(int version) const { return query(0, size(), version); }
  Monoid query() const { return query(latest()); }
  template <std::predicate<Monoid> Pred>
  int bisect_right(int l, const Pred& pred, int version) const {
#if !defined(NDEBUG)
    assert(0 <= l && l <= size());
    assert(pred(id_()));
    assert(0 <= version && version <= latest());
#endif
    if (l < size()) [[likely]] {
      return bisectr(tree_[version], l, pred);
    } else {
      return size();
    }
  }
  template <std::predicate<Monoid> Pred>
  int bisect_right(int l, const Pred& pred) const {
    return bisect_right(l, pred, latest());
  }
  template <std::predicate<Monoid> Pred>
  int bisect_left(int r, const Pred& pred, int version) const {
#if !defined(NDEBUG)
    assert(0 <= r && r <= size());
    assert(pred(id_()));
    assert(0 <= version && version <= latest());
#endif
    if (r > 0) [[likely]] {
      return bisectl(tree_[version], r, pred);
    } else {
      return 0;
    }
  }
  template <std::predicate<Monoid> Pred>
  int bisect_left(int r, const Pred& pred) const {
    return bisect_left(r, pred, latest());
  }
  int latest() const { return static_cast<int>(tree_.size()) - 1; }
  int size() const { return size_; }
private:
  template <detail::input_iterator<Monoid> IIter>
  node_pointer build(IIter base, int n) const {
    node_pointer res = node_allocator::allocate(*base);
    if (n > 1) {
      res->left = build(base, n / 2);
      res->right = build(base + n / 2, (n + 1) / 2);
      update(res);
    }
    return res;
  }
  node_pointer build(int n, const Monoid& value) const {
    node_pointer res = node_allocator::allocate(value);
    if (n > 1) {
      res->left = build(n / 2, value);
      res->right = build((n + 1) / 2, value);
      update(res);
    }
    return res;
  }
  const_node_pointer search(const_node_pointer root, int pos) const {
    int l = 0, r = size();
    const_node_pointer cur = root;
    while (r - l > 1) {
      detail::assume(l + r >= 0);
      const int mid = (l + r) / 2;
      if (pos < mid) {
        cur = cur->left, r = mid;
      } else {
        cur = cur->right, l = mid;
      }
    }
    return cur;
  }
  template <typename G>
  node_pointer search_update(const_node_pointer root, int pos, G g) const {
    int l = 0, r = size();
    std::stack<node_pointer, std::vector<node_pointer>> stk;
    node_pointer tmp = node_allocator::clone(root);
    node_pointer cur = tmp;
    while (r - l > 1) {
      stk.push(cur);
      detail::assume(l + r >= 0);
      const int mid = (l + r) / 2;
      if (pos < mid) {
        cur->left = node_allocator::clone(cur->left);
        cur = cur->left, r = mid;
      } else {
        cur->right = node_allocator::clone(cur->right);
        cur = cur->right, l = mid;
      }
    }
    cur->value = g(cur->value);
    while (!stk.empty()) { update(stk.top()), stk.pop(); }
    return tmp;
  }
  Monoid product(const_node_pointer root, int tl, int tr, int l, int r) const {
    if (l <= tl && tr <= r) { return root->value; }
    detail::assume(tl + tr >= 0);
    const int mid = (tl + tr) / 2;
    Monoid resl = l < mid ? product(root->left, tl, mid, l, r) : id_();
    Monoid resr = mid < r ? product(root->right, mid, tr, l, r) : id_();
    return binop_(std::move(resl), std::move(resr));
  }
  template <typename Pred>
  int bisectl(const_node_pointer root, int r, const Pred& pred) const {
    Monoid now = id_();
    auto fn = [&](auto fn, const_node_pointer cur, int tl, int tr) -> int {
      if (tr <= r) {
        if (Monoid tmp = binop_(cur->value, now); pred(tmp)) {
          now = std::move(tmp);
          return tl;
        }
        if (tr - tl == 1) { return tr; }
      }
      detail::assume(tl + tr >= 0);
      const int mid = (tl + tr) / 2;
      const int resr = mid < r ? fn(fn, cur->right, mid, tr) : mid;
      return resr != mid ? resr : fn(fn, cur->left, tl, mid);
    };
    return fn(fn, root, 0, size());
  }
  template <typename Pred>
  int bisectr(const_node_pointer root, int l, const Pred& pred) const {
    Monoid now = id_();
    auto fn = [&](auto fn, const_node_pointer cur, int tl, int tr) -> int {
      if (l <= tl) {
        if (Monoid tmp = binop_(now, cur->value); pred(tmp)) {
          now = std::move(tmp);
          return tr;
        }
        if (tr - tl == 1) { return tl; }
      }
      detail::assume(tl + tr >= 0);
      const int mid = (tl + tr) / 2;
      const int resl = l < mid ? fn(fn, cur->left, tl, mid) : mid;
      return resl != mid ? resl : fn(fn, cur->right, mid, tr);
    };
    return fn(fn, root, 0, size());
  }
  void update(node_pointer root) const {
    root->value = binop_(root->left->value, root->right->value);
  }
  void add_latest(node_pointer new_root) { tree_.emplace_back(new_root); }
  int size_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
  std::vector<node_pointer> tree_;
};
}  // namespace algo