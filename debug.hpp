#pragma once

#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
// makecode
#if defined(_DEBUG)
#  define PRINT(format, ...)                       \
    do {                                           \
      fmt::print("[{}:{}] ", __FILE__, __LINE__);  \
      fmt::print(FMT_STRING(format), __VA_ARGS__); \
      fmt::print("\n");                            \
    } while (0)
#  define DEBUG(x) \
    do { fmt::print("[{}:{}] {} = {}\n", __FILE__, __LINE__, #x, x); } while (0)
#else
#  define PRINT(...) (void)0
#  define DEBUG(x)   (void)0
#endif