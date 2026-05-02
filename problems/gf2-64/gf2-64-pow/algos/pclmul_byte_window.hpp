#pragma once
// PCLMUL + 4-bit window + Frobenius byte table。
//
// pow(a, e) は a と e 両方が変数なので linear_map で完全 byte 化はできない。
// しかし「multi-step squaring (Frobenius)」部分は GF(2)-線型なので byte table
// 化できる。これにより:
//   - 64 個の sq をシーケンシャルに繋ぐ部分が消える (64 × 2 cycle = 128 cycle)
//   - 代わりに 15 個の frob_4 (= a → a^{16}) を byte table で適用
//     (各 8 lookups, ~10 cycles, 15 回の chain で 150 cycles)
//
// 加えて 4-bit window で precompute (a^0..a^15) しておけば、e の各 nibble を
// 1 mul で消化できる。
//
// 期待: pclmul_pdep より速い、特に arm-clang (PDEP fallback が遅い) で大きく勝つ。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_RUN [[gnu::target("pclmul")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pow_byte_window {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// Frobenius 4-step (a → a^{2^4} = a^{16}) byte table
inline u64 FROB4_BYTE[8][256];
inline bool inited = false;

PCLMUL_FN void init_frob4_table() {
 if (inited) return;
 inited = true;
 u64 col[64];
 for (int j = 0; j < 64; ++j) {
  u64 v = u64(1) << j;
  for (int k = 0; k < 4; ++k) v = sq(v);  // 4 sqs = ^{2^4}
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

[[gnu::always_inline]] inline u64 frob4(u64 a) {
 return FROB4_BYTE[0][u8(a)]       ^ FROB4_BYTE[1][u8(a >>  8)]
      ^ FROB4_BYTE[2][u8(a >> 16)] ^ FROB4_BYTE[3][u8(a >> 24)]
      ^ FROB4_BYTE[4][u8(a >> 32)] ^ FROB4_BYTE[5][u8(a >> 40)]
      ^ FROB4_BYTE[6][u8(a >> 48)] ^ FROB4_BYTE[7][u8(a >> 56)];
}

PCLMUL_FN u64 pow(u64 a, u64 e) {
 if (e == 0) return 1;
 // Precompute T[i] = a^i for i = 0..15 (14 muls)
 u64 T[16];
 T[0] = 1;
 T[1] = a;
 #pragma GCC unroll 14
 for (int i = 2; i < 16; ++i) T[i] = mul(T[i-1], a);

 // Find top nonzero nibble of e
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

} // namespace gf2_64_pow_byte_window

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window::pow;
  using gf2_64_pow_byte_window::init_frob4_table;
  init_frob4_table();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
