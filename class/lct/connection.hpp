#pragma once

#include "../../internal/base/typing.hpp"
#include "../../internal/lct.hpp"
/*
@internal/lct.hpp
@internal/base/typing.hpp
*/

// makecode
namespace algo::link_cut {
namespace impl {
struct connectivity_node {
  int children[2]{0};
  int parent{0};
  bool flip{0};
};
}  // namespace impl
class connectivity_handler
    : public detail::basic_lct<connectivity_handler, impl::connectivity_node> {
  using self_type = connectivity_handler;
  using node_type = impl::connectivity_node;
  using base_type = detail::basic_lct<self_type, node_type>;
  friend class detail::basic_lct<self_type, node_type>;
public:
  connectivity_handler() noexcept = default;
  explicit connectivity_handler(int n) noexcept : base_type{n} {}
private:
  /**
   * \brief Pushing down the lazy tag to the children of u
   *
   * \param u
   */
  void push_down(int u) {
    if (tree_[u].flip) {
      std::swap(tree_[u].children[0], tree_[u].children[1]);
      const int l = tree_[u].children[0];
      const int r = tree_[u].children[1];
      tree_[u].flip ^= true;
      tree_[l].flip ^= true;
      tree_[r].flip ^= true;
    }
  }
};
}  // namespace algo::link_cut
