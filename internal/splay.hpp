#pragma once
#include <cstdint>
#include <memory>
#include <memory_resource>

#include "./base/def.hpp"
#include "./base/typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo::detail {
/**
 * \brief Basic splay tree node implementaion. This isn't meant to be used as
 * such. You must specialize: \brief 1. push_down(u) -> push lazy aggregate
 * information to the children. (u is not null) \brief 2. push_up(u) ->
 * recompute the aggregation.(u is not null) \brief 3. Node should contain the
 * member variable parent, children and size.
 * \brief 4. make(key) for creating node instance.
 * \brief 5. node.key for key type
 *
 * \tparam Node
 */
template <typename Derived, typename Node>
class basic_splay_tree {
public:
  using node_pointer = Node *;
  using const_node_pointer = const Node *;
  using key_type = typename Node::key_type;
  using pool_resource_type = std::pmr::monotonic_buffer_resource;
  basic_splay_tree() noexcept = default;
  // Declare move constructor to make the derived move-only type.
  // Which also delete the copy constructor and copy assignment operator.
  basic_splay_tree(basic_splay_tree &&tree) noexcept
      : root_{std::exchange(tree.root_, null())} {}
  // Move assignment operator
  basic_splay_tree &operator=(basic_splay_tree &&rhs) noexcept {
    root_ = std::exchange(rhs.root_, null());
    return *this;
  }
  /**
   * \brief Reroot the root to given node.
   *
   * \param u != null()
   */
  void reroot(node_pointer u) {
#if !defined(NDEBUG)
    assert(u != null());
#endif
    splay(u);
    root_ = u;
  }
  /**
   * \brief Get the address of node at position i. Which can be useful for
   * retrieving the position of node later.
   *
   * \param i 0 <= i <= size()
   * \return
   */
  node_pointer address(int i) {
#if !defined(NDEBUG)
    assert(0 <= i && i <= size());
#endif
    if (i == size()) { return null(); }
    const auto node = find_by_index(i);
    reroot(node);
    return node;
  }
  /**
   * \brief Get the index of node u.
   *
   * \param u u must be null() or valid address.
   * \return
   */
  int index(node_pointer u) {
    if (u == null()) { return size(); }
    reroot(u);
    return size(u->children[0]);
  }
  /**
   * \brief Get key of the node with index i.
   *
   * \param i
   * \return
   */
  key_type get(int i) {
#if !defined(NDEBUG)
    assert(0 <= i && i < size());
#endif
    return get(find_by_index(i));
  }
  /**
   * \brief Get key of the node with address u.
   *
   * \param u != null() is required
   * \return
   */
  key_type get(node_pointer u) {
#if !defined(NDEBUG)
    assert(u != null());
#endif
    reroot(u);
    return u->self;
  }
  /**
   * \brief Set new_key to the node with index i. Equivalent to set(address(i),
   * new_key)
   *
   * \param i
   * \param new_key
   */
  void set(int i, key_type new_key) {
#if !defined(NDEBUG)
    assert(0 <= i && i < size());
#endif
    set(find_by_index(i), std::move(new_key));
  }
  /**
   * \brief Set the new key to the node with address u.
   *
   * \param u
   * \param new_key
   */
  void set(node_pointer u, key_type new_key) {
#if !defined(NDEBUG)
    assert(u != null());
#endif
    reroot(u);
    u->self = std::move(new_key);
    update(u);
  }
  /**
   * \brief Insert the new node with given key immediately before the node with
   * index i. If i == size(), push the node at the end of the sequence.
   *
   * \param i
   * \param key
   * \return
   */
  node_pointer insert(int i, key_type key) {
#if !defined(NDEBUG)
    assert(0 <= i && i <= size());
#endif
    return insert(i, allocate(std::move(key)));
  }
  /**
   * \brief Insert the new node with given key immediately before the node with
   * address u. If u == null(), push the node at the end of the sequence.
   *
   * \param u
   * \param key
   * \return
   */
  node_pointer insert(node_pointer u, key_type key) {
    return insert(u, allocate(std::move(key)));
  }
  /**
   * \brief Insert new_node immediately before the index i and then return the
   * address of it. Equivalent to insert(address(i), new_node)
   *
   * \param i
   * \param new_node
   */
  node_pointer insert(int i, node_pointer new_node) {
#if !defined(NDEBUG)
    assert(0 <= i && i <= size());
#endif
    // To ensure the find_by_index(i) is not a bug.
    if (i < size()) {
      const auto node = find_by_index(i);
      // Make it the root.
      splay(node);
      insert_before(node, new_node);
      // Don't forget to splay it.
      reroot(new_node);
      return new_node;
    }
    // i == size() case
    return push_back(new_node);
  }
  /**
   * \brief Insert new_node immediately before the element of address u, and
   * then return the address of new_node. Equivalent to insert(index(u),
   * new_node)
   *
   * \param u u must be valid or null(). If u is null(), insert at the end.
   * \param new_node != nullptr
   * \return
   */
  node_pointer insert(node_pointer u, node_pointer new_node) {
#if !defined(NDEBUG)
    assert(new_node != null());
#endif
    if (u != null()) {
      // Make it the root.
      splay(u);
      insert_before(u, new_node);
      // Don't forget to splay it.
      reroot(new_node);
      return new_node;
    }
    return push_back(new_node);
  }
  /**
   * \brief Erase the node at position i. The remaining elements will be
   * concatenated.
   *
   * \param i 0 <= i < size() should hold.
   */
  void erase(int i) {
#if !defined(NDEBUG)
    assert(0 <= i && i < size());
#endif
    // Since 0 <= i < size(), find_by_index(i) is not a bug.
    erase(find_by_index(i));
  }
  /**
   * \brief Erase the node u. The remaining elements will be concatenated.
   *
   * \param u must be valid node of the tree.
   */
  void erase(node_pointer u) {
#if !defined(NDEBUG)
    assert(u != null());
#endif
    // Make node the root of the tree
    splay(u);
    const auto l = detach_discard(u, 0);
    const auto r = detach_discard(u, 1);
    root_ = concat(l, r);
  }
  /**
   * \brief Erase the elements lies in [l,r).
   *
   * \param l 0 <= l <= r <= size() should hold.
   * \param r
   */
  void erase(int l, int r) {
#if !defined(NDEBUG)
    assert(0 <= l && l <= r && r <= size());
#endif
    // To ensure that access_range(l,r) is not a bug.
    if (l < r) {
      // node != nullptr
      const auto node = access_range(l, r);
      const auto nxt = node->parent;
      // This check is required. Because access_range(l,r) might return the
      // root_ which has no parent after all.
      if (nxt != null()) {
        detach(nxt, side(node));
        splay(nxt);
      }
      // Set root_ to the valid root of the tree (or null)
      root_ = nxt;
    }
  }
  /**
   * \brief Extract node at index i and return the address of it.
   *
   * \param i
   * \return
   */
  node_pointer extract(int i) {
#if !defined(NDEBUG)
    assert(0 <= i && i < size());
#endif
    // Since i < size(), find_by_index(i) can be directly called
    return extract(find_by_index(i));
  }
  /**
   * \brief Extract node u from the tree and return it.
   *
   * \param u
   * \return
   */
  node_pointer extract(node_pointer u) {
#if !defined(NDEBUG)
    assert(u != null());
#endif
    splay(u);
    const auto l = detach(u, 0);
    const auto r = detach(u, 1);
    root_ = concat(l, r);
    return u;
  }
  /**
   * \brief Create the new node with given key and push it to the end of the
   * sequence.
   *
   * \param key
   * \return
   */
  node_pointer push_back(key_type key) {
    return push_back(allocate(std::move(key)));
  }
  /**
   * \brief Push new_node to the null of the sequence and return the address of
   * it.
   *
   * \param new_node != nullptr
   */
  node_pointer push_back(node_pointer new_node) {
    if (!empty()) [[likely]] {
      const auto node = root_;
      // nxt points to the last element of the sequence.
      const auto nxt = right_most(node);
      // nxt->children[1] == nullptr
      attach(nxt, new_node, 1);
      splay(new_node);
    }
    root_ = new_node;
    return new_node;
  }
  /**
   * \brief Create the new node with given key and prepend it onto the begging
   * of the sequence.
   *
   * \param key
   * \return
   */
  node_pointer push_front(key_type key) {
    return push_front(allocate(std::move(key)));
  }
  /**
   * \brief Prepend new_node to the beginning of the sequence and return the
   * address of it.
   *
   * \param new_node != null()
   */
  node_pointer push_front(node_pointer new_node) {
    if (!empty()) [[likely]] {
      const auto node = root_;
      // nxt points to the first element of the sequence.
      const auto nxt = left_most(node);
      // nxt->children[0] == nullptr
      attach(nxt, new_node, 0);
      splay(new_node);
    }
    root_ = new_node;
    return new_node;
  }
  /**
   * \brief Pop the last element of the sequence.
   *
   */
  void pop_back() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    const auto node = root_;
    // This points to the last element of the sequence
    const auto last = right_most(node);
    // Make it the root
    splay(last);
    // Detach left child of it.
    const auto nxt = detach_discard(last, 0);
    // nxt is already the root so no need to splay it.
    root_ = nxt;
  }
  /**
   * \brief Pop the first element of the sequence
   *
   */
  void pop_front() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    const auto node = root_;
    // This points to the first element of the sequence
    const auto first = left_most(node);
    // Make it the root
    splay(first);
    // Detach right child of it
    const auto nxt = detach_discard(first, 1);
    // nxt is already the root so no need to splay it.
    root_ = nxt;
  }
  /**
   * \brief Return the first node less than or equal to the given key assuming
   * the tree is ordered set. Otherwise the behavior is undefined.
   *
   * \tparam Compare
   * \param key
   * \param comp
   * \return
   */
  template <
      std::strict_weak_order<key_type, key_type> Compare = std::less<key_type>>
  node_pointer lower_bound(const key_type &key, const Compare &comp = {}) {
    node_pointer res = null();
    node_pointer last = root_;
    node_pointer cur = root_;
    while (cur != null()) {
      push(cur);
      last = cur;
      if (!comp(cur->self, key)) {
        // go left
        res = cur;
        cur = cur->children[0];
      } else {
        cur = cur->children[1];
      }
    }
    if (last != null()) { reroot(last); }
    return res;
  }
  template <
      std::strict_weak_order<key_type, key_type> Compare = std::less<key_type>>
  node_pointer upper_bound(const key_type &key, const Compare &comp = {}) {
    node_pointer res = null();
    node_pointer last = root_;
    node_pointer cur = root_;
    while (cur != null()) {
      push(cur);
      last = cur;
      if (comp(key, cur->self)) {
        res = cur;
        cur = cur->children[0];
      } else {
        cur = cur->children[1];
      }
    }
    if (last != null()) { reroot(last); }
    return res;
  }
  void clear() { root_ = null(); }
  int size() const { return size(root_); }
  bool empty() const { return size() == 0; }
  node_pointer null() const { return nullptr; }
protected:
  // // Extract elements lies in [l,r).
  // // if l == r, nullptr will be returned.
  // node_pointer extract(int l, int r) {
  //   // To ensure that access_range(l,r) is not a bug.
  //   if (l < r) [[likely]] {
  //     // node != nullptr
  //     const auto node = access_range(l, r);
  //     const auto nxt = node->parent;
  //     // This check is required. Because access_range(l,r) might return the
  //     // root_ which has no parent after all.
  //     if (nxt != nullptr) {
  //       detach(nxt, side(node));
  //       splay(nxt);
  //     }
  //     // Set root_ to the valid root of the tree (or null)
  //     root_ = nxt;
  //     return node;
  //   }
  //   // when the interval is empty
  //   return nullptr;
  // }
  /**
   * \brief Return the node which is responsible for the non-empty range [l,r).
   * After the method returns the previous root is not valid anymore. You should
   * call splay against the returned node to ensure the
   * subsequent calls of find_by_index(i) is not a bug.
   *
   * \param l 0 <= l < r <= size() should hold
   * \param r
   * \return Node which is responsible for the range. This node is not splayed.
   * The node is either root_, or the node which will become root by just only
   * rotation.
   */
  node_pointer access_range(int l, int r) const {
    const auto left = l != 0 ? find_by_index(l - 1) : null();
    const auto right = r != size() ? find_by_index(r) : null();
    if (left != null() && right != null()) {
      splay(left);
      splay(right);
      if (left->parent != right) { rotate(left); }
      // now left->parent == right
      return left->children[1];
    } else if (left != null()) {
      splay(left);
      return left->children[1];
    } else if (right != null()) {
      splay(right);
      return right->children[0];
    }
    // When l == 0 && r == size().
    return root_;
  }
  /**
   * \brief Find the node at index i. This method assumes that 'this->root_'
   * points to the valid root of the tree before this method is called.
   *
   * \param i 0 <= i < size() must be satisfied.
   * \return Node at index i. The node is not splayed.
   */
  node_pointer find_by_index(int i) const {
    auto cur = root_;
    while (true) {
      push(cur);
      const int left = size(cur->children[0]);
      if (i == left) {
        return cur;
      } else if (i < left) {
        cur = cur->children[0];
      } else {
        i -= left + 1;
        cur = cur->children[1];
      }
    }
    unreachable();
  }
  /**
   * \brief Merge two trees with roots l and r. Both have to be the roots of
   * their trees or null. The previous root is not valid after the call.
   *
   * \param l
   * \param r
   * \return New root of the joined tree. The node is splayed.
   */
  node_pointer concat(node_pointer l, node_pointer r) const {
    if (l != null()) {
      const auto nxt = right_most(l);
      attach(nxt, r, 1);
      splay(nxt);
      return nxt;
    }
    return r;
  }
  /**
   * \brief Return the left most node from u.
   *
   * \param u != nullptr
   * \return Left most node. The node is not splayed.
   */
  node_pointer left_most(node_pointer u) const {
    push(u);
    while (u->children[0] != null()) {
      u = u->children[0];
      push(u);
    }
    return u;
  }
  /**
   * \brief Return the right most node from u.
   *
   * \param u != nullptr
   * \return Right most node. The node is not splayed.
   */
  node_pointer right_most(node_pointer u) const {
    push(u);
    while (u->children[1] != null()) {
      u = u->children[1];
      push(u);
    }
    return u;
  }
  /**
   * \brief Insert u just before the v in the BST order.
   * \brief The previous root is valid after the call.
   * \brief After the call, you should splay against u.
   *
   * \param v != nullptr
   * \param u
   */
  void insert_before(node_pointer v, node_pointer u) const {
    if (v->children[0] == null()) {
      attach(v, u, 0);
    } else {
      const auto nxt = right_most(v->children[0]);
      attach(nxt, u, 1);
    }
  }
  /**
   * \brief Insert u just after the v in the BST order.
   * \brief The previous root is valid after the call.
   * \brief After the call, you should splay against u.
   *
   * \param v != nullptr
   * \param u
   */
  void insert_after(node_pointer v, node_pointer u) const {
    if (v->children[1] == null()) {
      attach(v, u, 1);
    } else {
      const auto nxt = left_most(v->children[1]);
      attach(nxt, u, 0);
    }
  }
  // u must not be null
  int side(const_node_pointer u) const {
    const auto p = u->parent;
    return !is_root(u) ? p->children[1] == u : -1;
  }
  /**
   * \brief Check if u is the root.
   *
   * \param u != nullptr
   * \return
   * \return
   */
  bool is_root(const_node_pointer u) const { return u->parent == null(); }
  /**
   * \brief Push lazy values to its children.
   *
   * \param u
   */
  void push(node_pointer u) const { self()->push_down(u); }
  /**
   * \brief Recompute the aggregation value of u from its children.
   *
   * \param u != nullptr
   */
  void update(node_pointer u) const { self()->push_up(u); }
  /**
   * \brief Attach u to v's one of its children.
   *
   * \param v != nullptr
   * \param u
   * \param d
   */
  void attach(node_pointer v, node_pointer u, bool d) const {
    v->children[d] = u;
    if (u != null()) { u->parent = v; }
    update(v);
  }
  /**
   * \brief Detach one of v's children
   *
   * \param v != nullptr
   * \param d
   * \return
   */
  node_pointer detach(node_pointer v, bool d) const {
    const auto u = v->children[d];
    v->children[d] = null();
    if (u != null()) { u->parent = null(); }
    update(v);
    return u;
  }
  /**
   * \brief Detach v's one of its children. v's aggregate values will not
   * updated.
   *
   * \param v != nullptr
   * \param d
   * \return The node which is detached from u. This node is not splayed.
   */
  node_pointer detach_discard(node_pointer v, bool d) const {
    const auto u = v->children[d];
    if (u != null()) { u->parent = null(); }
    return u;
  }
  // assume already propagated
  void rotate(node_pointer u) const {
    const auto v = u->parent;
    const auto w = v->parent;
    const int du = side(u);
    const int dv = side(v);
    const auto b = u->children[!du];
    attach(v, b, du);
    attach(u, v, !du);
    if (w != null()) { attach(w, u, dv); }
    u->parent = w;
  }
  /**
   * \brief Make u the root of the tree.
   *
   * \param u != nullptr
   */
  [[gnu::hot]] void splay(node_pointer u) const {
    push(u);
    while (!is_root(u) && !is_root(u->parent)) {
      const auto v = u->parent;
      const auto w = v->parent;
      push(w);
      push(v);
      push(u);
      rotate(side(u) == side(v) ? v : u);
      rotate(u);
    }
    if (!is_root(u)) {
      push(u->parent);
      push(u);
      rotate(u);
    }
  }
  /**
   * \brief Return the number of nodes in the subtree of node or 0 if node is
   * null.
   *
   * \param node
   * \return
   */
  int size(const_node_pointer node) const {
    return node != null() ? node->size : 0;
  }
  /**
   * \brief Allocate node based on the key.
   *
   * \param key
   * \return
   */
  node_pointer allocate(key_type key) const noexcept {
    void *mem = pool_.allocate(sizeof(Node), alignof(Node));
    return ::new (mem) Node{self()->make(std::move(key))};
  }
  auto self() const { return static_cast<const Derived *>(this); }
  auto self() { return static_cast<Derived *>(this); }
  /**
   * \brief Default push_down implementation. Which is no-op
   *
   */
  void push_down(node_pointer) const {}  // no-op
  /**
   * \brief Default push_up implementation. Which is no-op
   *
   */
  void push_up(node_pointer) const {}  // no-op
  node_pointer root_{nullptr};
  static inline pool_resource_type pool_;
};
}  // namespace algo::detail
