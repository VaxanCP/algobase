#pragma once

#include <array>
#include <cassert>
#include <stack>
#include <vector>

#include "../internal/base/typing.hpp"

/*
@internal/base/typing.hpp
*/

// makecode
namespace algo {
template <std::semiregular Monoid,
          detail::function<Monoid(Monoid, Monoid)> BinOp,
          detail::function<Monoid(void)> Id>
class swag_queue {
  using element_type = std::array<Monoid, 2>;
  using stack_type = std::stack<element_type, std::vector<element_type>>;
public:
  swag_queue() noexcept = default;
  swag_queue(BinOp binop, Id id) noexcept
      : front_{}, back_{}, binop_{std::move(binop)}, id_{std::move(id)} {}
  Monoid query() const {
    Monoid prod_f = !front_.empty() ? front_.top()[1] : id_();
    Monoid prod_b = !back_.empty() ? back_.top()[1] : id_();
    return binop_(std::move(prod_f), std::move(prod_b));
  }
  void push(Monoid value) {
    Monoid prod = !back_.empty() ? binop_(back_.top()[1], value) : value;
    back_.push({std::move(value), std::move(prod)});
  }
  void pop() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    // If front stack is empty, push all elements in back stack into the front
    // stack.
    if (front_.empty()) {
      front_.push({back_.top()[0], back_.top()[0]});
      back_.pop();
      while (!back_.empty()) {
        Monoid v = back_.top()[0];
        front_.push({v, binop_(v, front_.top()[1])});
        back_.pop();
      }
    }
    front_.pop();
  }
  bool empty() const { return front_.empty() && back_.empty(); }
  int size() const { return static_cast<int>(front_.size() + back_.size()); }
private:
  stack_type front_;
  stack_type back_;
  [[no_unique_address]] BinOp binop_;
  [[no_unique_address]] Id id_;
};
}  // namespace algo