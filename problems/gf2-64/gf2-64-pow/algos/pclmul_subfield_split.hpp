#pragma once
// Subfield decomposition pow:
//
//   F_{2^64}^* ≅ F_{2^16}^* × Q (norm-1 subgroup)
//   α = β · γ where β ∈ F_{2^16}^*, γ ∈ Q
//   α^e = β^{e mod 65535} · γ^{e mod M}    (M = (2^64-1)/65535 = 2^48 + 2^32 + 2^16 + 1)
//
// β^e は F_{2^16} log/exp で爆速 (~30 cycles)
// γ^e は 48-bit pow (元の 64-bit より ~25% 短い)
// ただし decompose の overhead (~50 cycles) と pow 削減 (~40 cycles) が拮抗するため
// 実機での効果は微妙な可能性が高い。実測で検証。
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

namespace gf2_64_pow_subfield_split {
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
inline u64 FROB4_BYTE[8][256];
inline u16 PW16[65536], LN16[65536];

// M^{-1} mod 65535 を事前計算 (compile-time):
//   M = 2^48 + 2^32 + 2^16 + 1 mod 65535
//   2^16 ≡ 1 mod 65535 なので M ≡ 1+1+1+1 = 4 mod 65535
//   ⇒ 4 · u ≡ 1 mod 65535、 u = 16384 (= (65535+1)/4)
constexpr u32 M_INV_MOD_65535 = 16384;
inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // frob16 byte table (16 sqs)
 {
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < 16; ++k) v = sq(v);
   col[j] = v;
  }
  for (int p = 0; p < 8; ++p) {
   for (int b = 0; b < 256; ++b) {
    u64 v = 0;
    for (int bit = 0; bit < 8; ++bit) {
     if ((b >> bit) & 1) v ^= col[p * 8 + bit];
    }
    FROB16_BYTE[p][b] = v;
   }
  }
 }
 // frob4 byte table for byte_window pow on γ
 {
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < 4; ++k) v = sq(v);
   col[j] = v;
  }
  for (int p = 0; p < 8; ++p) {
   for (int b = 0; b < 256; ++b) {
    u64 v = 0;
    for (int bit = 0; bit < 8; ++bit) {
     if ((b >> bit) & 1) v ^= col[p * 8 + bit];
    }
    FROB4_BYTE[p][b] = v;
   }
  }
 }
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

PCLMUL_FN u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}
PCLMUL_FN u64 frob4(u64 a) {
 return FROB4_BYTE[0][u8(a)]       ^ FROB4_BYTE[1][u8(a >>  8)]
      ^ FROB4_BYTE[2][u8(a >> 16)] ^ FROB4_BYTE[3][u8(a >> 24)]
      ^ FROB4_BYTE[4][u8(a >> 32)] ^ FROB4_BYTE[5][u8(a >> 40)]
      ^ FROB4_BYTE[6][u8(a >> 48)] ^ FROB4_BYTE[7][u8(a >> 56)];
}

// F_{2^16} subfield 元 (poly basis 64-bit) → 16-bit
PCLMUL_FN u16 extract_f16(u64 N_poly) {
 const u64 N_nim = gf2_64_basis::poly_to_nim(N_poly);
 return u16(N_nim);
}
// 16-bit F_{2^16} 元 → 64-bit poly basis
PCLMUL_FN u64 embed_f16(u16 n) {
 return gf2_64_basis::nim_to_poly(u64(n));
}
// F_{2^16} 内の log/exp ベースの pow (in poly basis)
PCLMUL_FN u64 f16_pow(u64 b_poly, u64 exp_value) {
 const u16 b16 = extract_f16(b_poly);
 if (b16 == 0) return 0;
 const u32 e_red = u32(exp_value % 65535);
 if (e_red == 0) return embed_f16(1);
 const u32 log_b = LN16[b16];
 const u32 log_result = u32((u64(log_b) * e_red) % 65535);
 return embed_f16(PW16[log_result]);
}

// γ ∈ F_{2^64} に対する byte_window pow
// 注: e_high = e mod M can be up to ~2^49 (M = 2^48 + 2^32 + 2^16 + 1 > 2^48)、
//     なので 12 nibbles では足りない。安全側で 16 nibbles (64-bit) フルレンジ。
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
 // Step 1: N(α) = α · α^{2^16} · α^{2^32} · α^{2^48}
 const u64 a16 = frob16(a);
 const u64 a32 = frob16(a16);
 const u64 a48 = frob16(a32);
 const u64 N = mul(mul(a, a16), mul(a32, a48));
 // Step 2: β = N^{u} where u = M^{-1} mod 65535 = 16384
 //   β は subfield 元なので 16-bit 表現でログ計算
 const u16 N16 = extract_f16(N);
 if (N16 == 0) {
  // α が norm 0、すなわち α 自身が 0 (普通は起こらないが)
  return 0;
 }
 const u32 log_N = LN16[N16];
 const u32 log_beta = u32((u64(log_N) * M_INV_MOD_65535) % 65535);
 const u16 beta16 = PW16[log_beta];
 const u64 beta_poly = embed_f16(beta16);
 // Step 3: γ = α · β^{-1}, β^{-1} は log/exp で
 const u32 log_beta_inv = (65535u - log_beta) % 65535u;
 const u16 beta_inv16 = PW16[log_beta_inv];
 const u64 beta_inv_poly = embed_f16(beta_inv16);
 const u64 gamma = mul(a, beta_inv_poly);
 // Step 4: e を分解 (e = e_low + 65535 * e_high の形ではなく、CRT で:
 //          β^e = β^{e mod 65535}, γ^e = γ^{e mod M})
 //   M = (2^64-1)/65535、 mod M は 64-bit 演算 (e は 64-bit)
 constexpr u64 M_VAL = (~u64(0)) / 65535u;  // = (2^64-1)/65535
 const u32 e_low  = u32(e % 65535);
 const u64 e_high = e % M_VAL;
 // Step 5: β^{e_low}, γ^{e_high}
 const u64 beta_pow = (e_low == 0) ? 1ull
   : embed_f16(PW16[u32((u64(log_beta) * e_low) % 65535)]);
 const u64 gamma_pow = pow_byte_window(gamma, e_high);
 // Step 6: 結合
 return mul(beta_pow, gamma_pow);
}

} // namespace gf2_64_pow_subfield_split

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split::pow;
  using gf2_64_pow_subfield_split::init_tables;
  init_tables();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
