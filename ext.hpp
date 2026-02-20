#pragma once
#include <functional>

#include "./internal/base/typing.hpp"
/*
@internal/base/typing.hpp
*/
// makecode
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

namespace algo {
namespace pb_ds = __gnu_pbds;
template <std::semiregular Key,
          std::relation<Key, Key> Compare = std::less<Key>>
using ordered_set =
    pb_ds::tree<Key, pb_ds::null_type, Compare, pb_ds::rb_tree_tag,
                pb_ds::tree_order_statistics_node_update>;
template <std::semiregular Key, std::semiregular Value,
          std::relation<Key, Key> Compare = std::less<Key>>
using ordered_map = pb_ds::tree<Key, Value, Compare, pb_ds::rb_tree_tag,
                                pb_ds::tree_order_statistics_node_update>;
template <std::semiregular Key, std::semiregular Value,
          detail::function<size_t(Key)> Hash = std::hash<Key>,
          std::equivalence_relation<Key, Key> Eq = std::equal_to<Key>>
using hash_map = pb_ds::gp_hash_table<Key, Value, Hash>;
template <std::semiregular Key,
          detail::function<size_t(Key)> Hash = std::hash<Key>,
          std::equivalence_relation<Key, Key> Eq = std::equal_to<Key>>
using hash_set = pb_ds::gp_hash_table<Key, pb_ds::null_type, Hash, Eq>;
}  // namespace algo