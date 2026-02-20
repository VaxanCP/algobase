#pragma once

#include <cassert>
#include <vector>

#include "./base/def.hpp"
#include "./base/typing.hpp"

/*
@internal/base/typing.hpp
@internal/base/def.hpp
*/

// makecode
namespace algo::detail {
/**
 * \brief Base class for link cut tree implementation. To customize behavior,
 * you must specialize(All nodes are assumed to be non-null unless explicitly
 * stated otherwise):
 * 1. push_down(u) for pushing down the lazy value of u to his left and right
 * children.
 * 2. push_up(u) for recomputing aggregation values from u's left and right
 * child.
 * 3. fix_link(v, u) for fixing subtree aggregation values and lazy invariants
 * when linking u to v's right child (u was previously not a virtual children of
 * v).
 * 4. fix_cut(v, u) for fixing subtree aggregation values and lazy invariants
 * when cutting u from v. (u is previously v's left child)
 * 5. fix_rotate(y, x, b) for fixing subtree aggregation values and lazy
 * invariants when rotating around x.(b maybe null)
 * 6. fix_remove(v, c) for fixing subtree aggregation values and lazy invariants
 * when removing c from v's right child and then make it v's one of virtual
 * children.(c maybe null)
 * 7. fix_attach(v, u) for fixing subtree aggregation values and lazy invariants
 * when attaching u to v's right child. (u maybe null or u was previously one of
 * virtual children of v).
 * 8. implement splay node struct. It must have children[2], parent, and flip
 * member variables to ensure that operation can be performed successfully.
 * \tparam Derived
 * \tparam Node Node aggregation.
 */
template <typename Derived, typename Node>
class basic_lct {
public:
  basic_lct() noexcept = default;
  explicit basic_lct(int n, const Node& def = Node{}) noexcept
      : tree_(n + 1, def) {}
  /**
   * \brief Make u as the root of the tree u belongs to.
   *
   * \param u A vertex
   */
  void reroot(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    do_reroot(u + 1);
  }
  /**
   * \brief Connect v and u. After the operation, the previous parent-child
   * relation ship is no longer valid (i.e, the root of the tree changes)
   *
   * \param v
   * \param u
   */
  void link(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
    assert(!connected(v, u));
#endif
    do_link(v + 1, u + 1);
  }
  /**
   * \brief Disconnecet u and v.  After the operation, the previous
   * parent-child relation ship is no longer valid (i.e, the root of the
   * tree changes)
   *
   * \param v
   * \param u
   */
  void cut(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    do_cut(v + 1, u + 1);
  }
  /**
   * \brief Get the parent of u if exists. Otherwise return Null.
   *
   * \param u
   * \return int
   */
  int parent(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return get_parent(u + 1) - 1;
  }
  /**
   * \brief Get the root of the tree to which u belongs
   *
   * \param u
   * \return int
   */
  int root(int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
#endif
    return get_root(u + 1) - 1;
  }
  /**
   * \brief Return lca of u and v or null if u and v is not connected
   *
   * \param v
   * \param u
   * \return int
   */
  int lca(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    return get_lca(v + 1, u + 1) - 1;
  }
  /**
   * \brief Check if u and v are connected
   *
   * \param v
   * \param u
   * \return true
   * \return false
   */
  bool connected(int v, int u) {
#if !defined(NDEBUG)
    assert(0 <= u && u < num_nodes());
    assert(0 <= v && v < num_nodes());
#endif
    return reachable(v + 1, u + 1);
  }
  /**
   * \brief Get the number of nodes in the tree.
   *
   * \return int
   */
  int num_nodes() const { return static_cast<int>(tree_.size()) - 1; }
  /**
   * \brief Return the indicator for node with null address.
   *
   * \return int
   */
  int null() const { return Null - 1; }
protected:
  /**
   * \brief Perform rerooting.
   *
   * \param u
   */
  void do_reroot(int u) {
    access(u);
    tree_[u].flip ^= true;
    push(u);
  }
  /**
   * \brief Pass an edge between v and u assuming they are previously not
   * connected. \brief v and u must not be connected otherwise it is no longer a
   * forest.
   *
   * \param v
   * \param u
   */
  void do_link(int v, int u) {
    do_reroot(u);
    access(v);
    tree_[v].children[1] = u;
    tree_[u].parent = v;
    self()->fix_link(v, u);
    update(v);
  }
  void do_cut(int v, int u) {
    do_reroot(u);
    access(v);
    if (tree_[v].children[0] == u) {
      tree_[v].children[0] = Null;
      tree_[u].parent = Null;
      self()->fix_cut(v, u);
      update(v);
    }
  }
  int get_parent(int u) {
    access(u);
    return tree_[u].children[0];
  }
  int get_root(int u) {
    access(u);
    while (tree_[u].children[0] != Null) {
      u = tree_[u].children[0];
      push(u);
    }
    splay(u);
    return u;
  }
  int get_lca(int v, int u) {
    if (u == v) { return u; }
    access(v);
    int w = Null, z = u;
    do {
      splay(z);
      toggle(z, w);
      w = z, z = tree_[z].parent;
    } while (z != Null);
    splay(u);
    return tree_[v].parent != Null ? w : Null;
  }
  bool reachable(int v, int u) { return get_lca(v, u) != Null; }
  // Make u as the root of its represented tree.
  void access(int u) {
    splay(u);
    toggle(u, Null);
    while (tree_[u].parent != Null) {
      const int w = tree_[u].parent;
      splay(w);
      toggle(w, u);
      rotate(u);
    }
  }
  int side(int u) const {
    const int p = tree_[u].parent;
    if (tree_[p].children[0] == u) { return 0; }
    if (tree_[p].children[1] == u) { return 1; }
    return -1;
  }
  // Check if u is the root of aux tree
  bool is_root(int u) const { return side(u) == -1; }
  // Switch preferred child
  void toggle(int v, int u) {
    const int c = tree_[v].children[1];
    self()->fix_remove(v, c);
    self()->fix_attach(v, u);
    tree_[v].children[1] = u;
    tree_[u].parent = v;
    update(v);
  }
  void rotate(int u) {
    const int v = tree_[u].parent;
    const int w = tree_[v].parent;
    const int du = side(u);
    const int dv = side(v);
    const int b = tree_[u].children[!du];
    // Fix invariants
    self()->fix_rotate(v, u, b);
    attach(v, b, du);
    attach(u, v, !du);
    if (dv != -1) { attach(w, u, dv); }
    tree_[u].parent = w;
  }
  // Attach u as a v's one of children (depending on d)
  void attach(int v, int u, bool d) {
    tree_[v].children[d] = u;
    tree_[u].parent = v;
    update(v);
  }
  // u must not be null
  void update(int u) { self()->push_up(u); }
  // u must not be null
  void push(int u) { self()->push_down(u); }
  // Splay u to be the root of aux tree
  void splay(int u) {
    push(u);
    while (!is_root(u) && !is_root(tree_[u].parent)) {
      const int v = tree_[u].parent;
      const int w = tree_[v].parent;
      push(w);
      push(v);
      push(u);
      rotate(side(v) == side(u) ? v : u);
      rotate(u);
    }
    if (!is_root(u)) {
      push(tree_[u].parent);
      push(u);
      rotate(u);
    }
  }
  /**
   * \brief No-op. You should specialize the behavior of the method such that
   * rerooting will be successfully performed.
   *
   * \param u
   */
  void push_down(int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of push_up in your derived class.
   *
   */
  void push_up(int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of fix_link in your derived class.
   *
   */
  void fix_link(int, int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of fix_cut in your derived class.
   *
   */
  void fix_cut(int, int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of fix_rotate in your derived class.
   *
   */
  void fix_rotate(int, int, int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of fix_remove in your derived class.
   *
   */
  void fix_remove(int, int) {}
  /**
   * \brief No-op. This will be called instead if you don't specialize the
   * behavior of fix_attach in your derived class.
   *
   */
  void fix_attach(int, int) {}
  /**
   * \brief Cast 'this' to the derived class pointer.
   *
   * \return
   */
  const Derived* self() const { return static_cast<const Derived*>(this); }
  /**
   * \brief Cast 'this' to the derived class pointer.
   *
   * \return
   */
  Derived* self() { return static_cast<Derived*>(this); }
  static constexpr int Null = 0;  // Alias for the null node.
  std::vector<Node> tree_;
};
}  // namespace algo::detail