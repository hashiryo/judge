#pragma once
// 多項式 GF(2)[x] / P(x) での拡張 Euclid 互除法による逆元計算 (P = x^64 + x^4 + x^3 + x + 1)。
//
// 不変量: s·a ≡ u (mod P), t·a ≡ v (mod P)
// 初期値: u = P (deg 64), s = 0; v = a (deg < 64), t = 1
// 各反復: 高次側を低次側で「シフト+XOR」で削る (deg(u) - deg(v) bit シフト)
//        対応する s, t も同じシフト+XOR で更新
//        最終的に v = gcd(a, P) = 1 → t ≡ a^{-1} (mod P)
//
// u, v は最大 deg 64 ⇒ u128 で保持。
// s, t は最大 deg ≈ 64 (deg s + deg r ≈ deg P) ⇒ u128、最後に mod P 還元。
//
// 反復回数は ~64 回 (毎反復で degree が確実に下がる)。
// 1 反復のコストは clz + シフト + XOR ≈ 数十 cycle 程度。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"  // mul は使わないが既存 algo と同じヘッダ集合に揃える
namespace gf2_64_extgcd {
using u128= __uint128_t;
inline int clz128(u128 x) {
 u64 hi= u64(x >> 64);
 return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll(u64(x));
}
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// u128 (deg ≤ 127 多項式) を mod P (= x^64 + 0x1B) で u64 に還元。
inline u64 reduce_modP(u128 x) {
 u64 h= x >> 64, d= h ^ (h << 1);
 return u64(x) ^ RED[h >> 60] ^ d ^ (d << 3);
}
inline u64 inv(u64 a) {
 if(!a) return 0;

 u128 s= 0, t= 1;
 u64 u= 0x1Bull ^ (a << (__builtin_clzll(a) + 1));  // P-a*x^shift
 s^= t << (__builtin_clzll(a) + 1);
 u64 v= a;
 while(v != 1) {
  int du= 127 - clz128(u);
  int dv= 63 - __builtin_clzll(v);
  if(du < dv) {
   u64 tmp= u;
   u= v;
   v= tmp;
   tmp= s;
   s= t;
   t= tmp;
   int td= du;
   du= dv;
   dv= td;
  }
  int shift= du - dv;
  u^= v << shift;
  s^= t << shift;
 }
 return reduce_modP(t);
}
inline u64 mul_naive(u64 a, u64 b) {
 // div = mul(a, inv(b)) のための単純 mul。本ファイル内で他の依存を増やさないために
 // bit-by-bit clmul + reduce で実装 (reference.hpp と同じ)。
 u64 lo= 0, hi= 0;
 for(int i= 0; i < 64; ++i) {
  if((b >> i) & 1) {
   lo^= a << i;
   if(i) hi^= a >> (64 - i);
  }
 }
 return reduce_modP((u128)hi << 64 | lo);
}
}  // namespace gf2_64_extgcd
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= gf2_64_extgcd::mul_naive(as[i], gf2_64_extgcd::inv(bs[i]));
  return ans;
 }
};
