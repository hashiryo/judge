#pragma once
// =============================================================================
// GF(2^64) ≅ GF(2)[x] / (x^64 + x^4 + x^3 + x + 1) の基本演算 (PCLMUL 経路)。
// 「**現状ベストの mul / sq / pow / reduce**」を提供する。
//
// === 運用ポリシー ===
// このファイルは **dynamic** で、 より速い実装が見つかったら更新される。
//
// 利用側のルール: 「自分の比較対象になる演算」は本ファイルから取り出してはいけない。
//
//   gf2-64-mul/algos/*    ← mul / sq を本ファイルから取らない (自己完結)
//   gf2-64-sq/algos/*     ← sq を本ファイルから取らない    (自己完結)
//   gf2-64-div/algos/*    ← inv を本ファイルから取らない   (mul / sq は OK)
//   gf2-64-pow/algos/*    ← pow を本ファイルから取らない   (mul / sq は OK)
//   gf2-64-sqrt/algos/*   ← sqrt を本ファイルから取らない  (mul / sq は OK)
//   gf2-64-log/algos/*    ← (log は本ファイルに無いので無関係)
//
//   理由: 比較対象の演算が "current best" に slot-in されてしまうと、 baseline が
//   暗黙に強くなって意味のある比較ができなくなる。 各比較対象は algo 内に閉じる。
//
// 更新時の手順:
//   1. 当該 problem の algos/ で新バリアントを追加し比較
//   2. 速い側に確定したら本ファイルの該当関数を更新
//   3. 他の問題は building block として自動的に高速化される
//
// CE 対策メモ:
//   GCC の `#pragma GCC target("pclmul")` は clang では無視される (warning のみ)。
//   関数単位の `__attribute__((target("pclmul")))` で代替するが、`always_inline` で
//   呼ぶ親関数も同じ target を持っていないと
//   「target が違う関数を always_inline できない」と clang がエラーを出す。
//   → 連鎖する全関数 + 呼び出し側 (algos/pclmul.hpp) に target("pclmul") を付ける。
// =============================================================================
#include <array>
#include <cstdint>

#include "_common.hpp"
namespace gf2_64_pclmul {
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
PCLMUL inline u64 mul(u64 a, u64 b) {
 __m128i v= _mm_clmulepi64_si128(_mm_cvtsi64_si128(a), _mm_cvtsi64_si128(b), 0);
 u64 h= (u64)v[1], d= h ^ (h << 1);
 return (u64)v[0] ^ RED[h >> 60] ^ d ^ (d << 3);
}
inline u64 spread_bits(u32 a) {
#ifdef __BMI2__
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 // fallback: bit interleave for 32-bit input
 u64 x= a;
 x= (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x= (x | (x << 8)) & 0x00FF00FF00FF00FFull;
 x= (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
 x= (x | (x << 2)) & 0x3333333333333333ull;
 x= (x | (x << 1)) & 0x5555555555555555ull;
 return x;
#endif
}
inline u64 sq(u64 a) {
 u64 h= spread_bits(a >> 32), d= h ^ (h << 1);
 return spread_bits(a) ^ RED[h >> 60] ^ d ^ (d << 3);
}
// 累乗 (a^e) 二進展開。
[[gnu::target("pclmul")]] inline u64 pow(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) res= mul(res, a);
  a= sq(a);
  e>>= 1;
 }
 return res;
}
// Fermat-style inverse: a^(2^64 - 2)。より高速な Itoh-Tsujii は別 algo で試す。
[[gnu::target("pclmul")]] inline u64 inv(u64 a) { return pow(a, ~(u64)1); }
// sqrt: Frobenius (^2) の逆 = ^(2^63)。
[[gnu::target("pclmul")]] inline u64 sqrt(u64 a) { return pow(a, (u64)1 << 63); }
}  // namespace gf2_64_pclmul
