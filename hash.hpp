#pragma once

#include <chrono>
#include <cstdint>

#include "./internal/base/typing.hpp"

/*
@internal/base/typing.hpp
*/
// makecode
namespace algo {
class split_hash {
public:
  static inline const uint64_t fx =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  static uint64_t hash64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }
  static uint32_t hash32(uint32_t x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }
  size_t operator()(detail::qword_fittable auto x) const {
    return hash64(x ^ fx);
  }
  size_t operator()(detail::dword_fittable auto x) const {
    return hash32(x ^ static_cast<uint32_t>(fx));
  }
};
}  // namespace algo