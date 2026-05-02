#pragma once
// Norm-based inv に「Frobenius byte table を 3 段 (frob16/32/48) 持つ」並列化。
//
// pclmul_norm.hpp は frob16 byte table を 3 回チェイン適用するため依存チェーン上
// 24 lookups が直列。frob32 / frob48 の byte table も持つことで β1, β2, β3 を
// 独立に計算でき、CPU の OoO で並列実行されることを期待。
//
// メモリ: 3 × 16 KiB = 48 KiB byte tables (L1 ぎりぎり)。
// 計算量: 同じく 24 lookups 合計、ただし依存チェーン解消。
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

namespace gf2_64_pclmul_norm_parallel {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;

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

// 3 段の Frobenius byte tables: α → α^{2^16}, α^{2^32}, α^{2^48}
inline u64 FROB16_BYTE[8][256];
inline u64 FROB32_BYTE[8][256];
inline u64 FROB48_BYTE[8][256];
inline bool frob_inited = false;

PCLMUL_FN void init_frob_tables() {
 if (frob_inited) return;
 frob_inited = true;
 for (int stage = 0; stage < 3; ++stage) {
  const int reps = 16 * (stage + 1);
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < reps; ++k) v = sq(v);
   col[j] = v;
  }
  auto* tbl = (stage == 0) ? &FROB16_BYTE[0][0]
            : (stage == 1) ? &FROB32_BYTE[0][0]
                           : &FROB48_BYTE[0][0];
  for (int p = 0; p < 8; ++p) {
   for (int b = 0; b < 256; ++b) {
    u64 v = 0;
    for (int bit = 0; bit < 8; ++bit) {
     if ((b >> bit) & 1) v ^= col[p * 8 + bit];
    }
    tbl[p * 256 + b] = v;
   }
  }
 }
}

[[gnu::always_inline]] inline u64 apply_byte_table(const u64 tbl[8][256], u64 a) {
 return tbl[0][u8(a)]       ^ tbl[1][u8(a >>  8)]
      ^ tbl[2][u8(a >> 16)] ^ tbl[3][u8(a >> 24)]
      ^ tbl[4][u8(a >> 32)] ^ tbl[5][u8(a >> 40)]
      ^ tbl[6][u8(a >> 48)] ^ tbl[7][u8(a >> 56)];
}
[[gnu::always_inline]] inline u64 frob16(u64 a) { return apply_byte_table(FROB16_BYTE, a); }
[[gnu::always_inline]] inline u64 frob32(u64 a) { return apply_byte_table(FROB32_BYTE, a); }
[[gnu::always_inline]] inline u64 frob48(u64 a) { return apply_byte_table(FROB48_BYTE, a); }

PCLMUL_FN u64 frob_repeat_pdep(u64 a, int rep) {
 for (int i = 0; i < rep; ++i) a = sq(a);
 return a;
}

// F_{2^16} 内の inv: N^{-1} = N^{2^16 - 2}, addition chain on 15 = 8+4+2+1
PCLMUL_FN u64 inv_in_f16(u64 N) {
 const u64 T1 = N;
 const u64 T2 = mul(T1, sq(T1));
 const u64 T4 = mul(T2, frob_repeat_pdep(T2, 2));
 const u64 T8 = mul(T4, frob_repeat_pdep(T4, 4));
 u64 acc = mul(frob_repeat_pdep(T8, 4), T4);
 acc = mul(frob_repeat_pdep(acc, 2), T2);
 acc = mul(frob_repeat_pdep(acc, 1), T1);
 return sq(acc);
}

PCLMUL_FN u64 inv(u64 a) {
 // β1, β2, β3 を独立 byte table で並列計算 (依存なし → ILP 効くはず)
 const u64 b1 = frob16(a);
 const u64 b2 = frob32(a);
 const u64 b3 = frob48(a);
 const u64 beta = mul(mul(b1, b2), b3);
 const u64 N = mul(a, beta);
 return mul(beta, inv_in_f16(N));
}

} // namespace gf2_64_pclmul_norm_parallel

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_parallel::inv;
  using gf2_64_pclmul_norm_parallel::init_frob_tables;
  init_frob_tables();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = mul(as[i], inv(bs[i]));
  return ans;
 }
};
