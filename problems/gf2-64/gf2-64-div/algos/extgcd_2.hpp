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
#include "../../_shared/mul.hpp"  // mul は使わないが既存 algo と同じヘッダ集合に揃える

using gf2_64_pclmul::mul;
namespace gf2_64_extgcd {
inline u64 inv(u64 a) {
 if(!a) return 0;
 int shift= __builtin_clzll(a) + 1;
 u64 u= 0x1Bull ^ (a << shift);  // P-a*x^shift
 u64 t= 1, s= 1ull << shift;
 while(a != 1) {
  int du= 63 - __builtin_clzll(u);
  int dv= 63 - __builtin_clzll(a);
  if(du < dv) {
   u64 tmp= u;
   u= a;
   a= tmp;
   tmp= s;
   s= t;
   t= tmp;
   int td= du;
   du= dv;
   dv= td;
  }
  shift= du - dv;
  u^= a << shift;
  s^= t << shift;
 }
 return t;
}
}  // namespace gf2_64_extgcd
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], gf2_64_extgcd::inv(bs[i]));
  return ans;
 }
};
