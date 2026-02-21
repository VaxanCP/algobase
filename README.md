# 🚀 algobase

A Modern, Zero-Overhead, and Extremely Optimized C++20 Template Library for Competitive Programming.

## ⚠️ Disclaimer / Author's Note

Please note the following before using or referencing this library:

1. **Retired from CP:** I am no longer an active competitive programmer. 
2. **No Maintenance or Support:** I am not currently involved in the software industry. Therefore, I cannot respond to any issues, bug reports, pull requests, or questions regarding this code. This repository is provided "as is," and you must use it entirely at your own risk.

**【注意事項 / 免責事項】**

本ライブラリをご利用、または参照される前に以下の点をご承知おきください。

1. **競技プログラミングからの引退:** 筆者はすでに競技プログラミングの世界から引退しています。
2. **ノーサポート:** 現在はソフトウェア業界にも身を置いていないため、本コードに関する問題（Issue）の解決、バグ対応、質問への回答等のメンテナンスは一切行えません。完全に「現状渡し（As is）」であり、自己責任でのご利用をお願いいたします。

---

## 🔥 Philosophy

`algobase` は、競技プログラミングにおける「定数倍の暴力」と「モダンC++の美しさ」を両立させるためにゼロから設計されたライブラリです。

古き悪き `std::enable_if` による黒魔術（SFINAE）を完全に捨て去り、C++20 の `concept` と `requires` 節によって、安全で、美しく、そしてコンパイラに極限の最適化を許すアーキテクチャを構築しています。

## ✨ Key Features

* **SFINAE-Free Architecture:**
  長大で読解不能なコンパイルエラーとはお別れです。すべての型制約は `algo::detail::integer` などのクリーンなコンセプトによって記述され、設計者の意図を100%の純度でコードに反映しています。
* **Compile-Time Safety & Type Promotion:**
  乗算時のオーバーフローを防ぐ `imul_result` によるコンパイル時の自動型昇格や、意図しないバグの温床となる `bool` 型の数値演算からの無慈悲な排除など、致命的なミスをコンパイラレベルで弾き飛ばします。
* **Hardware-Level Optimization:**
  L1/L2キャッシュへのヒット率を劇的に向上させるメモリ圧縮、分岐予測をコントロールする `[[unlikely]]` の徹底、SIMD/AVX2を意識したデータ構造設計により、圧倒的な実行速度を誇ります。

## 💻 Quick Example

C++20のコンセプトにより、シグネチャはここまでクリーンになります。SFINAEの呪文は一切不要です。

```cpp
namespace algo {

// 整数型（boolを除く）のみを受け付ける安全で美しい設計
template <detail::integer Tp>
constexpr Tp totient(Tp n) {
#if !defined(NDEBUG)
  assert(n > 0);
#endif
  return detail::totient(n);
}

// 乗算オーバーフローをコンパイル時に完全に防ぐ型昇格アーキテクチャ
template <detail::signed_integer T1, detail::signed_integer T2>
constexpr auto euclid_algo(T1 a, T2 b)
    -> std::array<std::common_type_t<T1, T2>, 3> {
  using Tp = std::common_type_t<T1, T2>;
  return detail::euclid_algo<Tp>(a, b);
}

}  // namespace algo
