#pragma once
#include <cassert>
#include <limits>
#include <memory_resource>

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
 * \brief Container for set of non-negative integers.
 *
 * \tparam Integral
 * \tparam BitSize
 */
template <detail::integer Integral,
          int BitSize = std::numeric_limits<Integral>::digits>
class binary_trie {
  static_assert(0 < BitSize &&
                BitSize <= std::numeric_limits<Integral>::digits);
  struct trie_node {
    int size;
    trie_node* children[2];
  };
  struct trie_node_allocator {
    using node_type = trie_node;
    using node_pointer = trie_node*;
    using const_node_pointer = const trie_node*;
    using pool_resource_type = std::pmr::monotonic_buffer_resource;
    static node_pointer allocate() {
      void* mem = pool_.allocate(sizeof(node_type), alignof(node_type));
      return ::new (mem) node_type{};
    }
  private:
    static inline pool_resource_type pool_;
  };
  using node_allocator = trie_node_allocator;
  using node_type = typename node_allocator::node_type;
  using node_pointer = typename node_allocator::node_pointer;
  using const_node_pointer = typename node_allocator::const_node_pointer;
public:
  binary_trie() noexcept : root_{node_allocator::allocate()} {}
  /**
   * \brief Insert n nodes whose value is x and return the
   * number of node whose value is x
   *
   * \param x value to insert
   * \param n
   * \return int
   */
  int insert(Integral x, int n = 1) const {
    if (n == 0) { return 0; }
    node_pointer cur = root_;
    cur->size += n;
    for (int i = BitSize - 1; i >= 0; --i) {
      const bool b = (x >> i) & 1;
      if (cur->children[b] == nullptr) {
        cur->children[b] = node_allocator::allocate();
      }
      cur = cur->children[b];
      cur->size += n;
    }
    return cur->size;
  }
  /**
   * \brief Erase nodes whose value is x and return the
   * number of removed nodes
   *
   * \param x
   * \param n
   * \return int
   */
  int erase(Integral x, int n = 1) const {
    std::array<node_pointer, BitSize + 1> stk{};
    stk[BitSize] = root_;
    for (int i = BitSize - 1; i >= 0; --i) {
      const bool b = (x >> i) & 1;
      if (size(stk[i + 1]->children[b]) == 0) { return 0; }
      stk[i] = stk[i + 1]->children[b];
    }
    const int removed = std::min(stk[0]->size, n);
    for (int i = 0; i <= BitSize; ++i) { stk[i]->size -= removed; }
    return removed;
  }
  /**
   * \brief Erase all nodes with value x
   *
   * \param x
   * \return
   */
  int erase_all(Integral x) const {
    std::array<node_pointer, BitSize + 1> stk{};
    stk[BitSize] = root_;
    for (int i = BitSize - 1; i >= 0; --i) {
      const bool b = (x >> i) & 1;
      if (size(stk[i + 1]->children[b]) == 0) { return 0; }
      stk[i] = stk[i + 1]->children[b];
    }
    const int removed = stk[0]->size;
    for (int i = 0; i <= BitSize; ++i) { stk[i]->size -= removed; }
    return removed;
  }
  /**
   * \brief Return the number of value of which prefix
   * [l,BitSize) match the corresponding bits of x
   *
   * \param x
   * \return constexpr int
   */
  int count(Integral x, int l = 0) const {
    const_node_pointer cur = root_;
    for (int i = BitSize - 1; i >= l; --i) {
      const bool b = (x >> i) & 1;
      if (size(cur->children[b]) == 0) { return 0; }
      cur = cur->children[b];
    }
    return cur->size;
  }
  /**
   * \brief
   *
   * \param x
   * \return
   * \return
   */
  bool contains(Integral x) const { return count(x) != 0; }
  /**
   * \brief Return the kth maximum of x^y (zero indexed).
   * size > k and k >= 0 must be satisfied.
   *
   * \param x
   * \param k
   * \return constexpr Integral
   */
  Integral xor_max(Integral x, int k = 0) const {
#if !defined(NDEBUG)
    assert(size() > k && k >= 0);
#endif
    Integral res = 0;
    const_node_pointer cur = root_;
    int remain = k + 1;
    for (int i = BitSize - 1; i >= 0; --i) {
      const bool b = (x >> i) & 1;
      const int right = size(cur->children[!b]);
      if (right >= remain) {
        res |= Integral(1) << i;
        cur = cur->children[!b];
      } else {
        remain -= right;
        cur = cur->children[b];
      }
    }
    return res;
  }
  /**
   * \brief Return the kth minimum of x^y (zero indexed).
   * size > k and k >= 0 must be satisfied
   *
   * \param x
   * \param k
   * \return constexpr Integral
   */
  Integral xor_min(Integral x, int k = 0) const {
#if !defined(NDEBUG)
    assert(size() > k && k >= 0);
#endif
    Integral res = Integral(0);
    const_node_pointer cur = root_;
    int remain = k + 1;
    for (int i = BitSize - 1; i >= 0; --i) {
      bool b = (x >> i) & 1;
      const int left = size(cur->children[b]);
      if (left < remain) {
        res |= Integral(1) << i;
        remain -= left;
        b = !b;
      }
      cur = cur->children[b];
    }
    return res;
  }
  /**
   * \brief Return the number of nodes such that y ^ x < z
   * where y is the value of node
   *
   * \param x
   * \param z
   * \return constexpr int
   */
  int xor_count_lt(Integral x, Integral z) const {
    int res = 0;
    const_node_pointer cur = root_;
    for (int i = BitSize - 1; i >= 0; --i) {
      const bool b = (x >> i) & 1;
      const bool v = (z >> i) & 1;
      if (v != 0) { res += size(cur->children[b]); }
      if (size(cur->children[v ^ b]) == 0) { break; }
      cur = cur->children[v ^ b];
    }
    return res;
  }
  /**
   * \brief Return the number of nodes such that y ^ x <= z
   * where y is the value of node
   *
   * \param x
   * \param z
   * \return constexpr int
   */
  int xor_count_le(Integral x, Integral z) const {
    return xor_count_lt(x, z) + count(x ^ z);
  }
  /**
   * \brief Return the number of nodes such that y ^ x > z
   * where y is the value of node
   *
   * \param x
   * \param z
   * \return constexpr int
   */
  int xor_count_gt(Integral x, Integral z) const {
    return size() - xor_count_le(x, z);
  }
  /**
   * \brief Return the number of nodes such that y ^ x >= z
   * where y is the value of node
   *
   * \param x
   * \param z
   * \return constexpr int
   */
  int xor_count_ge(Integral x, Integral z) const {
    return size() - xor_count_lt(x, z);
  }
  /**
   * \brief Return the kth maximum value of all node in the
   * trie. size() > 0 must hold
   *
   * \return constexpr Integral
   */
  Integral max(int k = 0) const {
#if !defined(NDEBUG)
    assert(size() != 0);
#endif
    return xor_max(0, k);
  }
  /**
   * \brief Return the kth minimum value of all node in the
   * trie. size() > 0 must hold
   *
   * \return constexpr Integral
   */
  Integral min(int k = 0) const {
#if !defined(NDEBUG)
    assert(size() != 0);
#endif
    return xor_min(0, k);
  }
  /**
   * \brief Return the number of nodes whose value less than
   * x
   *
   * \param x
   * \return constexpr int
   */
  int count_lt(Integral x) const { return xor_count_lt(0, x); }
  /**
   * \brief Return the number of nodes whose value less than
   * or equal to x
   *
   * \param x
   * \return constexpr int
   */
  int count_le(Integral x) const { return xor_count_le(0, x); }
  /**
   * \brief Return the number of nodes whose value greater
   * than x
   *
   * \param x
   * \return constexpr int
   */
  int count_gt(Integral x) const { return xor_count_gt(0, x); }
  /**
   * \brief Return the number of nodes whose value greater
   * or equal to x
   *
   * \param x
   * \return constexpr int
   */
  int count_ge(Integral x) const { return xor_count_ge(0, x); }
  /**
   * \brief Return the number of nodes in the trie
   *
   * \return constexpr int
   */
  int size() const { return size(root_); }
  /**
   * \brief Check if trie is empty
   *
   * \return true
   * \return false
   */
  bool empty() const { return size() == 0; }
private:
  // node may be null
  int size(const_node_pointer node) const {
    return node != nullptr ? node->size : 0;
  }
  node_pointer root_;
};
}  // namespace algo
