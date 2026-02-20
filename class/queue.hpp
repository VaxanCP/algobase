#pragma once

#include <cassert>
#include <vector>

#include "../internal/base/typing.hpp"
/*
@internal/base/typing.hpp
*/

// makecode
namespace algo {
template <typename Tp>
class buffer_queue {
public:
  buffer_queue() noexcept : removed_{0}, buffer_{} {}
  void push(Tp v) { buffer_.push_back(std::move(v)); }
  void pop() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    removed_++;
  }
  const Tp& front() const {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    return buffer_[removed_];
  }
  Tp& front() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    return buffer_[removed_];
  }
  const Tp& back() const {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    return buffer_.back();
  }
  Tp& back() {
#if !defined(NDEBUG)
    assert(!empty());
#endif
    return buffer_.back();
  }
  void reserve(int n) {
#if !defined(NDEBUG)
    assert(n >= 0);
#endif
    buffer_.reserve(n);
  }
  bool empty() const { return size() == 0; }
  int size() const { return static_cast<int>(buffer_.size()) - removed_; }
private:
  int removed_;
  std::vector<Tp> buffer_;
};
}  // namespace algo