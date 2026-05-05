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

// constexpr で [2][256] テーブルを構築
constexpr auto EMBED_BYTE= []() {
 constexpr std::array<uint64_t, 16> SUBFIELD_BASIS= {0x0000000000000001ULL, 0x5fbfaec6aeac0002ULL, 0xb06c601895640004ULL, 0xb013b5277b7c0008ULL, 0xb5ebb915248a0010ULL, 0x109bb25b2c600020ULL, 0xbf3bd95bd4190040ULL, 0x0fc66342279b0080ULL, 0xb6418f5e57c50100ULL, 0xaa194bd4b83f0200ULL, 0x1b5217b4dcc70400ULL, 0xbb06fa73867a0800ULL, 0x006fd55b23331000ULL, 0x4ae8fb39198c2000ULL, 0xfbd141b29b4f4000ULL, 0x1d9ce1776be78000ULL};
 std::array<std::array<uint64_t, 256>, 2> T{};
 for(int half= 0; half < 2; ++half) {
  for(int i= 0; i < 256; ++i) {
   uint64_t v= 0;
   for(int b= 0; b < 8; ++b) {
    if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
   }
   T[half][i]= v;
  }
 }
 return T;
}();
inline u64 embed_idx(u16 idx) { return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)]; }
constexpr u64 M_VAL2= (u64(1) << 32) + 1;
constexpr u32 INV2_F16= 32768;  // = 2^{-1} mod 65535
constexpr u32 INV2_F17= 32768;  // = (-2)^{-1} mod 65537 (= 65535^{-1} mod 65537)
// G_b 生成元: h2 = 2^{(2^32+1)(2^16-1)} = 2^{0x0000ffff0000ffff} ∈ F_{2^32}, order 65537
constexpr u64 h2= 0x1c1e79669b95a7ce;

inline u16 LN_SIGMA[65536];
inline u16 PW_SIGMA_IDX[65536];

inline u64 PW_H2[65537];  // PW_H2[k] = h2^k ∈ F_{2^32} ⊂ F_{2^64}
// G_b の log table: counting-sort で N=65537 個の (key, log) を 2^17 バケットに振り分け、
// 各バケット内を線形走査。key の下位 17 bit はバケット index と一致するので保存不要 →
// vs 1 entry に upper47(key) | log17 を pack して 8 byte/entry に圧縮。
//   vs = 65537 × 8 = 512 KiB
//   os = 131073 × 4 ≈ 512 KiB
//   合計 1 MiB (open-addressing 版 1.5 MiB から ~33% 削減)
struct H2BucketLog {
 static constexpr u32 N= 65537;
 static constexpr u32 HASH_SIZE= 1u << 17;
 static constexpr u64 MASK= HASH_SIZE - 1;
 u32 os[HASH_SIZE + 1];  // bucket start (BSS zero 初期化)
 u64 vs[N];              // packed: (key & ~MASK) | log
 void build(u64 g) {
  // Pass 1: count → os[i+1]
  u64 cur= 1;
  for(u32 k= 0; k < N; ++k) {
   ++os[(cur & MASK) + 1];
   cur= mul(cur, g);
  }
  // Prefix sum: os[i] = start of bucket i, os[H_S] = N
  for(u32 i= 1; i <= HASH_SIZE; ++i) os[i]+= os[i - 1];
  // Pass 2: fill (os[bucket] を書き込みポインタとして使用、後で復元)
  cur= 1;
  for(u32 k= 0; k < N; ++k) {
   u64 logv= u64(k) * INV2_F17 % 65537u;
   u32 bucket= u32(cur & MASK);
   vs[os[bucket]++]= (cur & ~MASK) | logv;
   cur= mul(cur, g);
  }
  // Restore: fill 後 os[i] = old_os[i+1] になっているので 1 つ右シフト
  for(u32 i= HASH_SIZE; i > 0; --i) os[i]= os[i - 1];
  os[0]= 0;
 }
 u32 lookup(u64 y) const {
  u64 yh= y & ~MASK;
  u32 m= u32(y & MASK);
  for(u32 i= os[m]; i < os[m + 1]; ++i) {
   if((vs[i] & ~MASK) == yh) return u32(vs[i] & MASK);
  }
  return ~u32(0);
 }
};
inline H2BucketLog LN_H2;
u64 pow_byte_window(u64 g, u64 e) {
 u64 T[4]= {1, g};
 for(int i= 2; i < 4; ++i) T[i]= mul(T[i - 1], g);
 int top= 31 - (__builtin_clzll(e) >> 1);
 u64 acc= T[(e >> (2 * top)) & 3];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob2(acc);
  u8 chunk= (e >> (2 * i)) & 3;
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // F_{2^16} σ chain (低 16 bit 識別 — PEXT_MASK = 0xFFFF と同等)
 constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= u16(cur);
  LN_SIGMA[idx]= u32(k) * INV2_F16 % 65535u;  // log_σ(cur) · (65537)^{-1} を記録
  PW_SIGMA_IDX[k]= idx;
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
 PW_SIGMA_IDX[65535]= 1;
 // G_b chain: PW_H2[k] = h2^k for k = 0..65536
 cur= 1;
 for(u32 k= 0; k < 65537; ++k) {
  PW_H2[k]= cur;
  cur= mul(cur, h2);
 }
 // LN_H2 (counting-sort バケット) を別パスで構築
 LN_H2.build(h2);
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 const u32 q= u32(e / M_VAL2);
 if(!q) return pow_byte_window(a, e);
 const u64 r= e - u64(q) * M_VAL2;
 const u64 N2= mul(a, frob32(a));
 // x_a = N2(α)^{65537} = α · α^{2^16} · α^{2^32} · α^{2^48} ∈ G_a (= F_{2^16}^*)
 const u64 x_a= mul(N2, frob16(N2));
 // x_b = N2(α)^{65535} ∈ G_b (Itoh-Tsujii 4 mul + 4 frob)
 const u64 T3= mul(N2, sq(N2));          // N2^3
 const u64 T15= mul(T3, frob2(T3));      // N2^15
 const u64 T255= mul(T15, frob4(T15));   // N2^255
 const u64 x_b= mul(T255, frob8(T255));  // N2^65535
 const u16 L_a= LN_SIGMA[u16(x_a)];      // log_σ(x_a) · (65537)^{-1} = l_a mod 65535
 const u32 L_b= LN_H2.lookup(x_b);       // log_h2(x_b) · (65535)^{-1} = l_b mod 65537
 // l_a × q mod 65535, l_b × q mod 65537 を計算
 const u16 e_a= u64(q) * L_a % 65535u;
 const u32 e_b= u64(q) * L_b % 65537u;
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
