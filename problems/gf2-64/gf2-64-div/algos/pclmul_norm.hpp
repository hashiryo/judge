#pragma once
// Norm-based inversion (subfield reduction) on poly basis.
//
// 塔 F_{2^64} ⊃ F_{2^16} を活用。norm を使うと F_{2^64}-inv を F_{2^16}-inv に reduce できる:
//
//   N(α) = α · α^{2^16} · α^{2^32} · α^{2^48}  ∈ F_{2^16}  (Galois 共役の積は subfield)
//   α^{-1} = (α^{2^16} · α^{2^32} · α^{2^48}) · N(α)^{-1}
//
// 各 Frobenius α^{2^16} は GF(2)-線型なので 64×64 行列 → 8 byte tables (16 KiB) で
// 8 回 byte lookup に置換可能。16 PCLMUL-square より少ない命令数で済む。
//
// N は F_{2^16} subfield に住む。inv は F_{2^64} 上の Itoh-Tsujii で計算するが、
// N^{2^16 - 1} = 1 なので a^{-1} = a^{2^16 - 2} で済み、F_{2^15 - 1} までの addition chain
// で十分 (= 6 muls + 15 sqs)。
//
// 期待 ops: 10 muls + 15 sqs + 24 byte lookups
// (cf. pclmul_itoh_tsujii: 10 muls + 63 sqs)
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pclmul_norm {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;

// PDEP-based F_{2^64} squaring (= GF(2) では bit を 1 つおきに spread + reduce)
[[gnu::always_inline]] inline u64 spread_bits(u32 a) {
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

// Frobenius α → α^{2^16} の byte table (起動時に sq^16 で構築、8 byte tables × 256 entries)
inline u64 FROB16_BYTE[8][256];
inline bool frob_inited = false;

PCLMUL_FN void init_frob_table() {
 if (frob_inited) return;
 frob_inited = true;
 // basis vector x^j に Frobenius を 16 回適用 → α_j ∈ F_{2^64}
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

[[gnu::always_inline]] inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}

PCLMUL_FN u64 frob_repeat(u64 a, int rep) {
 for (int i = 0; i < rep; ++i) a = sq(a);
 return a;
}

// F_{2^16} 内の inv: N^{-1} = N^{2^16 - 2} = (N^{2^15 - 1})^2
// 15 = 8+4+2+1 で addition chain (Itoh-Tsujii 縮小版)
PCLMUL_FN u64 inv_in_f16(u64 N) {
 const u64 T1 = N;
 const u64 T2 = mul(T1, sq(T1));                  // a^{2^2 - 1}, 1 sq + 1 mul
 const u64 T4 = mul(T2, frob_repeat(T2, 2));      // a^{2^4 - 1}, 2 sqs + 1 mul
 const u64 T8 = mul(T4, frob_repeat(T4, 4));      // a^{2^8 - 1}, 4 sqs + 1 mul
 // T_15 = T_{8+4+2+1}: acc から 4,2,1 sqs ずつシフトして T_4, T_2, T_1 を掛ける
 u64 acc = mul(frob_repeat(T8, 4), T4);           // a^{2^12 - 1}, 4 sqs + 1 mul
 acc = mul(frob_repeat(acc, 2), T2);              // a^{2^14 - 1}, 2 sqs + 1 mul
 acc = mul(frob_repeat(acc, 1), T1);              // a^{2^15 - 1}, 1 sq + 1 mul
 return sq(acc);                                  // a^{2^16 - 2} = a^{-1}, 1 sq
}

PCLMUL_FN u64 inv(u64 a) {
 // β1, β2, β3 = α^{2^16}, α^{2^32}, α^{2^48}
 const u64 b1 = frob16(a);
 const u64 b2 = frob16(b1);
 const u64 b3 = frob16(b2);
 // β = b1 · b2 · b3
 const u64 beta = mul(mul(b1, b2), b3);
 // N = α · β  (∈ F_{2^16} ⊂ F_{2^64})
 const u64 N = mul(a, beta);
 // α^{-1} = β · N^{-1}
 return mul(beta, inv_in_f16(N));
}

} // namespace gf2_64_pclmul_norm

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm::inv;
  using gf2_64_pclmul_norm::init_frob_table;
  init_frob_table();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = mul(as[i], inv(bs[i]));
  return ans;
 }
};
