#pragma once

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

// makecode
namespace algo {
std::vector<int> prefix_function(std::string_view sv) {
  std::vector<int> pf(sv.length());
  for (int i = 1, j = 0; i < std::ssize(pf); ++i) {
    while (j != 0 && sv[i] != sv[j]) { j = pf[j - 1]; }
    if (sv[i] == sv[j]) { j++; }
    pf[i] = j;
  }
  return pf;
}
std::vector<int> z_function(std::string_view sv) {
  std::vector<int> zf(sv.length());
  for (int i = 1, l = 0, r = 0; i <= std::ssize(zf); i++) {
    if (i <= r) { zf[i] = std::min(zf[i - l], r - i + 1); }
    while (i + zf[i] < std::ssize(zf) && sv[zf[i]] == sv[i + zf[i]]) {
      zf[i]++;
    }
    if (i + zf[i] - 1 > r) {
      l = i;
      r = i + zf[i] - 1;
    }
  }
  return zf;
}
}  // namespace algo