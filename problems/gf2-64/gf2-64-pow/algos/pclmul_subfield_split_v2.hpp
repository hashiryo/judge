#pragma once
// pclmul_subfield_split.hpp の改良版:
//
// 改良 1: N(α) 計算で 3 frob16 シリアルチェーン (~30 cycle) を削除
//   - FROB16, FROB32, FROB48 の 3 つの byte table を持つ
//   - a16 = frob16(a), a32 = frob32(a), a48 = frob48(a) を全部 a から並列計算
//   - 依存解消で latency 短縮 (memory は 48 KiB、L1 borderline)
//
// 改良 2: tree 化された mul で N(α) の結合パスも短縮
//   - mul(mul(a, a16), mul(a32, a48)) は depth 2 tree で OK (元から)
//
// 改良 3: e mod 65535 を高速計算 (2^16 ≡ 1 なので桁和)
//   - e mod 65535 = (e0 + e1 + e2 + e3) mod 65535 where e_i = 16-bit chunks
//   - 64-bit % は ~10 cycle、桁和 + 1 reduce は ~5 cycle
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"
#include "../../_shared/basis_change.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pow_subfield_split_v2 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;
using u16 = unsigned short;
using u32 = unsigned;

PCLMUL_FN u64 spread_bits(u32 a) {
#if defined(__BMI2__) && !defined(USE_SIMDE)
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 u64 x = a;
 x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x = (x | (x <<  8)) & 0x00FF00FF00FF00FFull;
 x = (x | (x <<  4)) & 0x0F0F0F0F0F0F0F0Full;
 x = (x | (x <<  2)) & 0x3333333333333333ull;
 x = (x | (x <<  1)) & 0x5555555555555555ull;
 return x;
#endif
}
PCLMUL_FN u64 sq(u64 a) {
 const u64 lo = spread_bits(u32(a));
 const u64 hi = spread_bits(u32(a >> 32));
 return reduce(__m128i{(long long) lo, (long long) hi});
}

inline u64 FROB16_BYTE[8][256];
inline u64 FROB32_BYTE[8][256];
inline u64 FROB48_BYTE[8][256];
inline u64 FROB4_BYTE[8][256];
inline u16 PW16[65536], LN16[65536];

constexpr u32 M_INV_MOD_65535 = 16384;
inline bool inited = false;

PCLMUL_FN void build_frob_byte_table(int reps, u64 (&out)[8][256]) {
 u64 col[64];
 for (int j = 0; j < 64; ++j) {
  u64 v = u64(1) << j;
  for (int k = 0; k < reps; ++k) v = sq(v);
  col[j] = v;
 }
 for (int p = 0; p < 8; ++p) {
  for (int b = 0; b < 256; ++b) {
   u64 v = 0;
   for (int bit = 0; bit < 8; ++bit) {
    if ((b >> bit) & 1) v ^= col[p * 8 + bit];
   }
   out[p][b] = v;
  }
 }
}

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 build_frob_byte_table(4,  FROB4_BYTE);
 build_frob_byte_table(16, FROB16_BYTE);
 build_frob_byte_table(32, FROB32_BYTE);
 build_frob_byte_table(48, FROB48_BYTE);
 // F_{2^16} log/exp (Nimber 互換)
 PW16[0] = PW16[65535] = 1;
 for (int i = 1; i < 65535; ++i) {
  PW16[i] = u16((PW16[i-1] << 1) ^ (0x1681fu & u16(-(PW16[i-1] >= 0x8000u))));
 }
 constexpr u16 f2n[16] = {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u,
                          0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for (int i = 1; i < 65535; ++i) {
  u16 x = PW16[i], y = 0;
  for (; x; x &= x - 1) y ^= f2n[__builtin_ctz(x)];
  PW16[i] = y;
  LN16[y] = u16(i);
 }
 LN16[1] = 0;
}

PCLMUL_FN u64 apply_byte_table(const u64 (&t)[8][256], u64 a) {
 return t[0][u8(a)]       ^ t[1][u8(a >>  8)]
      ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)]
      ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)]
      ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
PCLMUL_FN u64 frob4 (u64 a) { return apply_byte_table(FROB4_BYTE,  a); }
PCLMUL_FN u64 frob16(u64 a) { return apply_byte_table(FROB16_BYTE, a); }
PCLMUL_FN u64 frob32(u64 a) { return apply_byte_table(FROB32_BYTE, a); }
PCLMUL_FN u64 frob48(u64 a) { return apply_byte_table(FROB48_BYTE, a); }

PCLMUL_FN u16 extract_f16(u64 N_poly) {
 return u16(gf2_64_basis::poly_to_nim(N_poly));
}
PCLMUL_FN u64 embed_f16(u16 n) {
 return gf2_64_basis::nim_to_poly(u64(n));
}

// e mod 65535 を桁和で高速計算 (2^16 ≡ 1 mod 65535 を利用)
[[gnu::always_inline]] inline u32 e_mod_65535(u64 e) {
 const u32 s = u32(e & 0xFFFF) + u32((e >> 16) & 0xFFFF)
             + u32((e >> 32) & 0xFFFF) + u32((e >> 48) & 0xFFFF);
 // s は最大 4 × 65535 = 262140 → mod 65535 で最大 2 回減算
 u32 r = s;
 if (r >= 65535) r -= 65535;
 if (r >= 65535) r -= 65535;
 if (r >= 65535) r -= 65535;
 return r;
}

PCLMUL_FN u64 pow_byte_window(u64 g, u64 e) {
 if (e == 0) return 1;
 u64 T[16];
 T[0] = 1; T[1] = g;
 #pragma GCC unroll 14
 for (int i = 2; i < 16; ++i) T[i] = mul(T[i-1], g);
 int top = 15;
 while (top > 0 && ((e >> (4 * top)) & 0xF) == 0) --top;
 u64 acc = T[(e >> (4 * top)) & 0xF];
 for (int i = top - 1; i >= 0; --i) {
  acc = frob4(acc);
  unsigned chunk = unsigned((e >> (4 * i)) & 0xF);
  if (chunk) acc = mul(acc, T[chunk]);
 }
 return acc;
}

PCLMUL_FN u64 pow(u64 a, u64 e) {
 if (e == 0) return 1;
 // ① N(α) 計算: 3 frob を独立 byte table で並列
 const u64 a16 = frob16(a);
 const u64 a32 = frob32(a);
 const u64 a48 = frob48(a);
 const u64 N = mul(mul(a, a16), mul(a32, a48));
 // ② β 抽出 (log/exp F_{2^16})
 const u16 N16 = extract_f16(N);
 if (N16 == 0) return 0;
 const u32 log_N = LN16[N16];
 const u32 log_beta = u32((u64(log_N) * M_INV_MOD_65535) % 65535);
 // ③ γ = α · β^{-1} via log/exp
 const u32 log_beta_inv = (65535u - log_beta) % 65535u;
 const u64 beta_inv_poly = embed_f16(PW16[log_beta_inv]);
 const u64 gamma = mul(a, beta_inv_poly);
 // ④ e の分解
 constexpr u64 M_VAL = (~u64(0)) / 65535u;
 const u32 e_low_red = e_mod_65535(e);
 const u64 e_high = e % M_VAL;
 // ⑤ β^{e_low}
 const u64 beta_pow = (e_low_red == 0) ? 1ull
   : embed_f16(PW16[u32((u64(log_beta) * e_low_red) % 65535)]);
 // ⑥ γ^{e_high}
 const u64 gamma_pow = pow_byte_window(gamma, e_high);
 // ⑦ 結合
 return mul(beta_pow, gamma_pow);
}

} // namespace gf2_64_pow_subfield_split_v2

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v2::pow;
  using gf2_64_pow_subfield_split_v2::init_tables;
  init_tables();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
