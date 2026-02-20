#pragma once

#include <concepts>
#include <functional>
#include <utility>
// makecode

namespace algo {
template <typename Fn>
class rec_lambda {
public:
  constexpr rec_lambda() noexcept = default;
  constexpr rec_lambda(Fn fn) noexcept : fn_{std::move(fn)} {}
  template <typename... Args>
  constexpr decltype(auto) operator()(Args&&... args) const {
    return fn_(*this, std::forward<Args>(args)...);
  }
  template <typename... Args>
  constexpr decltype(auto) operator()(Args&&... args) {
    return fn_(*this, std::forward<Args>(args)...);
  }
private:
  [[no_unique_address]] Fn fn_;
};
}  // namespace algo
