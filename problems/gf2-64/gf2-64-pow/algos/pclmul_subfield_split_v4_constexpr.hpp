#pragma once
// pclmul_subfield_split_v4_2 の preprocessing 完全 constexpr 化版:
//
// LN_SIGMA / PW_SIGMA_IDX (各 ~128 KiB) を natural 表現での σ chain で compile-time 構築。
// PCLMUL の constexpr mul は不要 — chain は 16-bit polynomial の (cur<<1) ^ (Q_LOW & -hi)
// で済むので 65535 周しても constexpr step 上限に余裕で収まる。
//
// 利点: init_tables 不要 → 初回 query の cold start 無し、コンパイラの定数畳み込みも狙える。
//
// 実装メモ: subfield_split_constexpr_2.hpp (div 側) と同じ手法。
//   - chain 生成元 BETA を nat 表現の "y" として使い 16-bit shift+XOR でループ。
//   - col[i] = u16(BETA^i) (poly basis) を hardcode、natural→low 変換 byte table を作る。
//   - SUBFIELD_BASIS (subfield 元の low → poly 64bit) は生成元非依存なので v4_2 と同じ。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v4_7 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob3;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// embed: low 16-bit subfield 識別子 → 64-bit poly 表現 (subfield 元の埋め込み)
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
// LN_SIGMA[u16(BETA^k)] = k, PW_SIGMA_IDX[k] = u16(BETA^k)
// natural chain で 65535 周しつつ T_lo/T_hi byte table で nat→low 変換し、両テーブル同時に埋める。
struct Tables {
 u16 LN_SIGMA[65536];
 u16 PW_SIGMA_IDX[65535];
};
constexpr auto TABLES= []() {
 // col[i] = u16(BETA^i) — BETA = 0x1f1af3ec55a22e02 の poly basis 累乗の低 16 bit
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
 // BETA の最小多項式は y^16 + y^5 + y^3 + y^2 + 1, lower bits = 0x002D
 u16 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u16 lo= T_lo[u8(cur)] ^ T_hi[cur >> 8];
  t.LN_SIGMA[lo]= u16(k);
  t.PW_SIGMA_IDX[k]= lo;
  cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 }
 t.LN_SIGMA[0]= 0;
 return t;
}();
u64 pow_byte_window(u64 g, u64 e) {
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
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
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
}  // namespace gf2_64_pow_subfield_split_v4_7
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v4_7::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
