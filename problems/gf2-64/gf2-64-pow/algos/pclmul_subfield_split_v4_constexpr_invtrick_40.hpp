#pragma once
// pclmul_subfield_split_v4_constexpr.hpp + popcount inversion trick:
//
//   popcount(e) > 32 のとき a^e = (a^{-1})^{~e} で計算 (~e = 2^64-1-e in u64)。
//   F_{2^64}^* は order 2^64-1 なので: (a^{-1})^{2^64-1-e} = a^{-(2^64-1-e)}
//                                    = a^{e + 1 - 2^64} ≡ a^e (mod 2^64-1) ✓
//   popcount(~e) = 64 - popcount(e) < 32 になるので pow_byte_window のコスト減。
//
// inv(a) は subfield_split_constexpr_2.hpp と同じく INV_LOW lookup 1 回で済ませる
// (LN/PW から導出する案もあるが lookup 2 回かかるので INV_LOW を持っておく)。
//
// 効果見込み: dense_e (popcount ~60) ケースで pow_byte_window の chunk 数が減る分速くなる。
//   トレードオフ: inv 計算 (~6 mul + 2 lookups) のオーバーヘッドが popcount 32 付近で
//   ブレイクイーブン。 popcount > 40 で明確にプラス。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v4_constexpr_invtrick {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob3;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u64 embed_idx(u16 idx) {
 static constexpr auto EMBED= []() {
  u64 SUBFIELD_BASIS[]= {0x0000000000000001ULL, 0x5fbfaec6aeac0002ULL, 0xb06c601895640004ULL, 0xb013b5277b7c0008ULL, 0xb5ebb915248a0010ULL, 0x109bb25b2c600020ULL, 0xbf3bd95bd4190040ULL, 0x0fc66342279b0080ULL, 0xb6418f5e57c50100ULL, 0xaa194bd4b83f0200ULL, 0x1b5217b4dcc70400ULL, 0xbb06fa73867a0800ULL, 0x006fd55b23331000ULL, 0x4ae8fb39198c2000ULL, 0xfbd141b29b4f4000ULL, 0x1d9ce1776be78000ULL};
  array<array<u64, 256>, 2> t{};
  for(int half= 0; half < 2; ++half)
   for(int i= 0; i < 256; ++i) {
    u64 v= 0;
    for(int b= 0; b < 8; ++b)
     if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
    t[half][i]= v;
   }
  return t;
 }();
 return EMBED[0][u8(idx)] ^ EMBED[1][idx >> 8];
}
struct Tables {
 u16 LN_SIGMA[65536];
 u16 PW_SIGMA_IDX[65535];
 u16 INV_LOW[65536];
};
constexpr auto TABLES= []() {
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {}, T_hi[256]= {};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j)
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
   }
  T_lo[v]= lo;
  T_hi[v]= hi;
 }
 Tables t{};
 // chain で LN_SIGMA / PW_SIGMA_IDX を埋める
 u16 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u16 lo= T_lo[u8(cur)] ^ T_hi[cur >> 8];
  t.LN_SIGMA[lo]= u16(k);
  t.PW_SIGMA_IDX[k]= lo;
  cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 }
 t.LN_SIGMA[0]= 0;
 // INV_LOW: σ^k inverse = σ^{65535-k} (順序 65535)
 //   PW_SIGMA_IDX[k] = u16(σ^k) なのでそれを使えば nat[] や T_lo/T_hi 再 lookup 不要。
 t.INV_LOW[t.PW_SIGMA_IDX[0]]= t.PW_SIGMA_IDX[0];  // σ^0 = 1 → INV[1] = 1
 for(u32 k= 1; k <= 32767; ++k) {
  u16 lo_k= t.PW_SIGMA_IDX[k];
  u16 lo_ik= t.PW_SIGMA_IDX[65535 - k];
  t.INV_LOW[lo_k]= lo_ik;
  t.INV_LOW[lo_ik]= lo_k;
 }
 return t;
}();
// a の F_{2^64}^* における逆元: subfield norm 経由 + INV_LOW lookup 1 回
inline u64 inv(u64 a) {
 u64 g= mul(frob16(a), mul(frob32(a), frob48(a)));
 u64 b= embed_idx(TABLES.INV_LOW[u16(mul(g, a))]);
 return mul(b, g);
}
inline u64 pow_byte_window(u64 g, u64 e) {
 u64 T[8]= {1, g};
 for(int i= 2; i < 8; ++i) T[i]= mul(T[i - 1], g);
 int top= 16 - ((__builtin_clzll(e) - 13) / 3);
 u64 acc= T[(e >> (3 * top)) & 7];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob3(acc);
  u8 chunk= (e >> (3 * i)) & 7;
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
inline u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 // popcount(e) > 32 のとき a^e = (a^{-1})^{~e}
 if(__builtin_popcountll(e) > 40) {
  a= inv(a);
  e= ~e;
  if(!e) return 1;  // e == 2^64-1 のケース → (a^{-1})^0 = 1
 }
 constexpr u64 M_VAL= (~u64(0)) / 65535u;
 const u16 q= e / M_VAL;
 if(!q) return pow_byte_window(a, e);
 const u64 N2= mul(a, frob32(a));
 const u16 N= mul(N2, frob16(N2));
 const u64 b= embed_idx(TABLES.PW_SIGMA_IDX[(u32(TABLES.LN_SIGMA[N]) * q) % 65535]);
 const u64 r= e - M_VAL * q;
 if(!r) return b;
 const u64 g= pow_byte_window(a, r);
 return mul(b, g);
}
}  // namespace gf2_64_pow_subfield_split_v4_constexpr_invtrick
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v4_constexpr_invtrick::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
