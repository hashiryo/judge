#pragma once
// 多項式 GF(2)[x] / P(x) での拡張 Euclid 互除法による逆元計算 (P = x^64 + x^4 + x^3 + x + 1)。
//
// 不変量: s·a ≡ u (mod P), t·a ≡ v (mod P)
// 初期値: u = P (deg 64), s = 0; v = a (deg < 64), t = 1
// 各反復: 高次側を低次側で「シフト+XOR」で削る (deg(u) - deg(v) bit シフト)
//        対応する s, t も同じシフト+XOR で更新
//        最終的に v = gcd(a, P) = 1 → t ≡ a^{-1} (mod P)
//
// u, v は最大 deg 64 ⇒ __uint128_t で保持。
// s, t は最大 deg ≈ 64 (deg s + deg r ≈ deg P) ⇒ __uint128_t、最後に mod P 還元。
//
// 反復回数は ~64 回 (毎反復で degree が確実に下がる)。
// 1 反復のコストは clz + シフト + XOR ≈ 数十 cycle 程度。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"  // mul は使わないが既存 algo と同じヘッダ集合に揃える
namespace gf2_64_extgcd {
[[gnu::always_inline]] inline int clz128(__uint128_t x) {
 u64 hi= u64(x >> 64);
 return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll(u64(x));
}
// __uint128_t (deg ≤ 127 多項式) を mod P (= x^64 + 0x1B) で u64 に還元。
[[gnu::always_inline]] inline u64 reduce_modP(__uint128_t x) {
 u64 lo= u64(x);
 u64 hi= u64(x >> 64);
 // hi · x^64 ≡ hi · 0x1B (mod P)。0x1B = 1 + x + x^3 + x^4 = bits {0,1,3,4}。
 __uint128_t prod= (__uint128_t)hi ^ ((__uint128_t)hi << 1) ^ ((__uint128_t)hi << 3) ^ ((__uint128_t)hi << 4);
 lo^= u64(prod);
 u64 hi2= u64(prod >> 64);  // ≤ 4 bit
 lo^= hi2 ^ (hi2 << 1) ^ (hi2 << 3) ^ (hi2 << 4);
 return lo;
}
inline u64 inv(u64 a) {
 if(!a) return 0;
 __uint128_t u= ((__uint128_t)1 << 64) | 0x1Bull;  // P
 __uint128_t v= a;
 __uint128_t s= 0, t= 1;
 while(v != 1) {
  int du= 127 - clz128(u);
  int dv= 127 - clz128(v);
  if(du < dv) {
   __uint128_t tmp= u;
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
 return reduce_modP((__uint128_t)hi << 64 | lo);
}
}  // namespace gf2_64_extgcd
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= gf2_64_extgcd::mul_naive(as[i], gf2_64_extgcd::inv(bs[i]));
  return ans;
 }
};
