#pragma once
// std headers
#include <compare>
#include <concepts>
#include <cstdint>
#include <iterator>
#include <type_traits>

// makecode
namespace algo::detail {
struct modnum_base {};
struct monostate {};
template <typename Tp>
concept integer = std::integral<Tp> && !(std::same_as<Tp, bool>);
template <typename Tp>
concept signed_integer = integer<Tp> && (Tp(-1) < Tp(0));
template <typename Tp>
concept unsigned_integer = integer<Tp> && (Tp(-1) > Tp(0));
template <typename Tp>
concept modular_integer = std::is_base_of_v<modnum_base, Tp>;
template <typename From, typename To>
concept sign_compatible_with =
    integer<From> && ((std::is_signed_v<From> && std::is_signed_v<To>) ||
                      (std::is_unsigned_v<From> && std::is_unsigned_v<To>));
template <typename Tp>
concept arithmetic = integer<Tp> || (std::is_floating_point_v<Tp>);
template <typename Tp>
concept qword_fittable = sizeof(Tp) <= 8UL;
template <typename Tp>
concept dword_fittable = qword_fittable<Tp> && sizeof(Tp) <= 4UL;
template <typename Tp>
concept word_fittable = dword_fittable<Tp> && sizeof(Tp) <= 2UL;
template <typename Tp>
concept byte_fittable = word_fittable<Tp> && sizeof(Tp) == 1UL;
template <typename Tp>
struct imul_result;
template <qword_fittable Tp>
  requires std::signed_integral<Tp>
struct imul_result<Tp> {
  using type = __int128_t;
};
template <qword_fittable Tp>
  requires std::unsigned_integral<Tp>
struct imul_result<Tp> {
  using type = __uint128_t;
};
template <dword_fittable Tp>
  requires std::signed_integral<Tp>
struct imul_result<Tp> {
  using type = int64_t;
};
template <dword_fittable Tp>
  requires std::unsigned_integral<Tp>
struct imul_result<Tp> {
  using type = uint64_t;
};
template <word_fittable Tp>
  requires std::signed_integral<Tp>
struct imul_result<Tp> {
  using type = int32_t;
};
template <word_fittable Tp>
  requires std::unsigned_integral<Tp>
struct imul_result<Tp> {
  using type = uint32_t;
};
template <byte_fittable Tp>
  requires std::signed_integral<Tp>
struct imul_result<Tp> {
  using type = int16_t;
};
template <byte_fittable Tp>
  requires std::unsigned_integral<Tp>
struct imul_result<Tp> {
  using type = uint16_t;
};
template <typename Tp>
using imul_result_t = typename imul_result<Tp>::type;
/// usage: function_traits<int(int,long)>
template <typename Sig>
struct function_traits;
template <typename Ret, typename... Args>
struct function_traits<Ret(Args...)> {
  using return_type = Ret;
  template <typename Fn>
  static constexpr bool same_as =
      std::is_invocable_r_v<return_type, Fn, Args...>;
};
template <typename Fn, typename Sig>
concept function = (function_traits<Sig>::template same_as<Fn>);
template <typename IIter, typename Value>
concept input_iterator = std::same_as<std::iter_value_t<IIter>, Value>;
template <typename OIter, typename Value>
concept output_iterator = std::output_iterator<OIter, Value>;
}  // namespace algo::detail
