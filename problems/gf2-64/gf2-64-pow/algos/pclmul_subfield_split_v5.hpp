#pragma once
// 境界 2^32 - 1 で α^e を分解する pow バリアント (v4_4 が境界 2^16 - 1 = 65535 だった発展):
//   α^e = N2(α)^q · α^r
//   N2(α) = α · α^{2^32} ∈ F_{2^32}^*
//   M_VAL2 = (2^64 - 1)/(2^32 - 1) = 2^32 + 1
//   q = e / M_VAL2 ∈ [0, 2^32 - 1], r = e mod M_VAL2 ∈ [0, 2^32]
//
// F_{2^32}^* の order = 2^32 - 1 = 65535 × 65537 を CRT 分解:
//   - G_a (order 65535) = F_{2^16}^*  → 既存 LN_SIGMA / PW_SIGMA_IDX で log/exp
//   - G_b (order 65537)               → 直接 hash log (pohlig_v7 風) + 配列 exp
//
// 抽出:
//   x_a = N2(α)^{65537} = α · α^{2^16} · α^{2^32} · α^{2^48} = N(α)  (v4_4 の N と同じ)
//   x_b = N2(α)^{65535}  (Itoh-Tsujii で 4 mul + 4 frob ≈ 60 cycle)
// CRT 補正係数 32768 = 2^{-1} mod 65535 = (-2)^{-1} mod 65537 (両方とも同じ値):
//   l_a = LN_SIGMA[x_a] × 32768 mod 65535   (∵ x_a = (α-part)^2)
//   l_b = LN_H2[x_b]    × 32768 mod 65537   (∵ x_b = (β-part)^{65535} = (β-part)^{-2})
// CRT 結合 32768 × 65537 + 32768 × 65535 = 2^32 ≡ 1 (mod 2^32 - 1) で α^e を再構成。
//
// メモリ:
//   既存 F_{2^16}: LN_SIGMA + PW_SIGMA_IDX = 256 KiB, EMBED_BYTE = 4 KiB
//   新   G_b:     PW_H2 = 65537 × 8 ≈ 512 KiB, LN_H2 hash (cap 131072) ≈ 1.5 MiB
//
// 利点: r ∈ [0, 2^32] なので byte_window pow が 8 chunks (v4_4 は 11 chunks)。
//       chunk_top_bit ケースなど r が大きい指数で効果。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v5 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::frob8;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
// σ^0..σ^15 (低 16 bit が線型独立、よって u16(cur) で識別可能 — v4_4 と同じ)
constexpr array<u64, 16> SIGMA_POW= {1ull, 11625825068197362226ull, 16726069499881173557ull, 1643712016162803251ull, 1139490258686435223ull, 2133510969792054338ull, 6902720057649047833ull, 6911235282440967732ull, 1609028687227055289ull, 17188600949917388119ull, 970962174813342580ull, 13487762679230420322ull, 17964985024964438593ull, 449778467375550729ull, 13547138401716955404ull, 5392541099451007413ull};
// EMBED_BYTE: 16-bit subfield idx → 64-bit poly basis (v4_4 と同じ Gauss-Jordan)
constexpr auto EMBED_BYTE= []() {
 array<array<u64, 256>, 2> embed_byte{};
 u32 M[16]= {}, Minv[16]= {};
 for(int r= 0; r < 16; ++r) {
  u32 row= 0;
  for(int i= 0; i < 16; ++i)
   if((SIGMA_POW[i] >> r) & 1) row|= u32(1) << i;
  M[r]= row;
  Minv[r]= u32(1) << r;
 }
 for(int col= 0; col < 16; ++col) {
  int pv= -1;
  for(int r= col; r < 16; ++r)
   if((M[r] >> col) & 1) {
    pv= r;
    break;
   }
  if(pv == -1) continue;
  if(pv != col) {
   u32 tm= M[col];
   M[col]= M[pv];
   M[pv]= tm;
   tm= Minv[col];
   Minv[col]= Minv[pv];
   Minv[pv]= tm;
  }
  for(int r= 0; r < 16; ++r)
   if(r != col && ((M[r] >> col) & 1)) {
    M[r]^= M[col];
    Minv[r]^= Minv[col];
   }
 }
 u64 contribution[16]= {};
 for(int k= 0; k < 16; ++k) {
  u64 v= 0;
  for(int i= 0; i < 16; ++i)
   if((Minv[i] >> k) & 1) v^= SIGMA_POW[i];
  contribution[k]= v;
 }
 for(int p= 0; p < 2; ++p)
  for(int b= 0; b < 256; ++b) {
   u64 v= 0;
   for(int bit= 0; bit < 8; ++bit)
    if((b >> bit) & 1) v^= contribution[p * 8 + bit];
   embed_byte[p][b]= v;
  }
 return embed_byte;
}();

constexpr u64 M_VAL2= (u64(1) << 32) + 1;
constexpr u32 INV2_F16= 32768;  // = 2^{-1} mod 65535
constexpr u32 INV2_F17= 32768;  // = (-2)^{-1} mod 65537 (= 65535^{-1} mod 65537)
// G_b 生成元: h2 = 2^{(2^32+1)(2^16-1)} = 2^{0x0000ffff0000ffff} ∈ F_{2^32}, order 65537
constexpr u64 H2_EXP= 0x0000ffff0000ffffull;

inline u16 LN_SIGMA[65536];
inline u16 PW_SIGMA_IDX[65536];

inline u64 PW_H2[65537];  // PW_H2[k] = h2^k ∈ F_{2^32} ⊂ F_{2^64}
struct H2Hash {
 static constexpr u64 CAP= 131072;  // 2 × 65537 を超える 2 冪 (load ~50%)
 static constexpr u64 MASK= CAP - 1;
 u64 keys[CAP];
 u32 vals[CAP];
 void clear() {
  for(u64 i= 0; i < CAP; ++i) keys[i]= ~u64(0);
 }
 void insert(u64 key, u32 val) {
  u64 h= (key * 0x9E3779B97F4A7C15ull) & MASK;
  while(keys[h] != ~u64(0)) h= (h + 1) & MASK;
  keys[h]= key;
  vals[h]= val;
 }
 u32 lookup(u64 key) const {
  u64 h= (key * 0x9E3779B97F4A7C15ull) & MASK;
  while(keys[h] != ~u64(0)) {
   if(keys[h] == key) return vals[h];
   h= (h + 1) & MASK;
  }
  return ~u32(0);  // not found (subfield 元が想定通りなら起きない)
 }
};
inline H2Hash LN_H2;

inline u64 embed_idx(u16 idx) { return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)]; }

u64 pow_byte_window(u64 g, u64 e) {
 if(!e) return 1;  // r==0 安全策
 u64 T[16]= {1, g};
#pragma GCC unroll 14
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], g);
 int top= 15 - (__builtin_clzll(e) >> 2);
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u32 chunk= u32((e >> (4 * i)) & 0xF);
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}

inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // F_{2^16} σ chain (低 16 bit 識別 — PEXT_MASK = 0xFFFF と同等)
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= u16(cur);
  LN_SIGMA[idx]= k;
  PW_SIGMA_IDX[k]= idx;
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
 PW_SIGMA_IDX[65535]= 1;
 // G_b chain: h2^k for k = 0..65536
 const u64 h2= pow_byte_window(2, H2_EXP);
 LN_H2.clear();
 cur= 1;
 for(u32 k= 0; k < 65537; ++k) {
  PW_H2[k]= cur;
  LN_H2.insert(cur, k);
  cur= mul(cur, h2);
 }
}

u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 const u32 q= u32(e / M_VAL2);
 if(!q) return pow_byte_window(a, e);
 const u64 r= e - u64(q) * M_VAL2;
 // x_a = N(α) = α · α^{2^16} · α^{2^32} · α^{2^48} ∈ G_a (= F_{2^16}^*)
 const u64 x_a= mul(mul(a, frob16(a)), mul(frob32(a), frob48(a)));
 // x_b = N2(α)^{65535} ∈ G_b (Itoh-Tsujii 4 mul + 4 frob)
 const u64 N2= mul(a, frob32(a));
 const u64 T3= mul(N2, sq(N2));               // N2^3
 const u64 T15= mul(T3, frob2(T3));           // N2^15
 const u64 T255= mul(T15, frob4(T15));        // N2^255
 const u64 x_b= mul(T255, frob8(T255));       // N2^65535
 const u32 L_a= LN_SIGMA[u16(x_a)];           // log_σ(x_a) = 2 · l_a mod 65535
 const u32 L_b= LN_H2.lookup(x_b);            // log_h2(x_b) = -2 · l_b mod 65537
 // l_a × q mod 65535, l_b × q mod 65537 を計算 (l_* = L_* × 32768 を吸収)
 const u32 q1= q % 65535u;
 const u32 q2= q % 65537u;
 const u32 e_a= u32((u64(L_a) * INV2_F16 % 65535u) * q1 % 65535u);
 const u32 e_b= u32((u64(L_b) * INV2_F17 % 65537u) * q2 % 65537u);
 const u64 b_a= embed_idx(PW_SIGMA_IDX[e_a]);  // SIGMA^{e_a} (= α-part^q)
 const u64 b_b= PW_H2[e_b];                    // h2^{e_b}    (= β-part^q)
 const u64 b= mul(b_a, b_b);
 if(!r) return b;
 const u64 g= pow_byte_window(a, r);
 return mul(b, g);
}
}  // namespace gf2_64_pow_subfield_split_v5
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v5::init_tables;
  using gf2_64_pow_subfield_split_v5::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
