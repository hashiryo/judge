#pragma once
// pclmul_norm_pext.hpp の memory 圧縮版:
// INV_PEXT を 64-bit poly (512 KiB) ではなく 16-bit idx (128 KiB) に変更し、
// idx → 64-bit poly の embed を 2 byte tables (4 KiB) で実現。
//
// embed は GF(2)-線型 → 16-bit input を low/high byte に分けて 2 lookup + XOR。
//
// メモリ: 元 512 KiB → 132 KiB (L2 余裕で fit)
// per inv: PEXT 1 + INV_PEXT_IDX 1 lookup + 2 EMBED byte lookup ≈ 18 cycle
//   (元 pext: PEXT 1 + INV_PEXT 1 lookup ≈ 12 cycle)
// memory pressure 軽減で実機性能はおそらく良くなる。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#define HAVE_PEXT 1
#else
#define PCLMUL_RUN
#define HAVE_PEXT 0
#endif

namespace gf2_64_pclmul_norm_pext_compact {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;
using gf2_64_pclmul::sq;
using u16 = unsigned short;
using u32 = unsigned;

constexpr u64 SIGMA = 0xa1573a4da2bc3a32ull;

inline u64 FROB16_BYTE[8][256];
inline u64 PEXT_MASK = 0;
inline u16 INV_PEXT_IDX[65536];   // PEXT idx → idx of inverse (16-bit, 128 KiB)
inline u64 EMBED_BYTE[2][256];    // 16-bit idx → 64-bit poly via 2 byte lookup (4 KiB)

#if !HAVE_PEXT
inline int PEXT_POS[16];
#endif
inline bool inited = false;

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
PCLMUL_FN u64 sq_pdep(u64 a) {
 const u64 lo = spread_bits(u32(a));
 const u64 hi = spread_bits(u32(a >> 32));
 return reduce(__m128i{(long long) lo, (long long) hi});
}

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // Frobenius byte table (16 sqs)
 {
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < 16; ++k) v = sq_pdep(v);
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
 // σ^0..σ^15 + Gauss 消去で線型独立な 16 bit
 u64 sigma_pow[16];
 sigma_pow[0] = 1;
 for (int i = 1; i < 16; ++i) sigma_pow[i] = mul(sigma_pow[i-1], SIGMA);
 u16 row_vec[64];
 for (int r = 0; r < 64; ++r) {
  u16 v = 0;
  for (int c = 0; c < 16; ++c) {
   if ((sigma_pow[c] >> r) & 1) v |= u16(1) << c;
  }
  row_vec[r] = v;
 }
 int picked[16]; int n_picked = 0;
 u16 basis[16] = {};
 for (int r = 0; r < 64 && n_picked < 16; ++r) {
  u16 v = row_vec[r];
  for (int k = 15; k >= 0 && v; --k) {
   if (!((v >> k) & 1)) continue;
   if (basis[k] == 0) {
    basis[k] = v;
    picked[n_picked++] = r;
    break;
   }
   v ^= basis[k];
  }
 }
 PEXT_MASK = 0;
 for (int i = 0; i < 16; ++i) PEXT_MASK |= (u64(1) << picked[i]);
#if !HAVE_PEXT
 for (int i = 0; i < 16; ++i) PEXT_POS[i] = picked[i];
#endif
 // σ^k chain で INV_PEXT_IDX を構築 + EMBED 用 contribution を取得
 std::vector<u64> PW_local(65535);
 std::vector<u32> PEXT_IDX_local(65535);
 PW_local[0] = 1;
 PEXT_IDX_local[0] = 0;
 for (int i = 0; i < 16; ++i) PEXT_IDX_local[0] |= 0;  // = 1 in F_{2^16} → idx of "1"
 // 実際 σ^0 = 1 の PEXT idx を計算
 {
  u32 idx0;
#if HAVE_PEXT
  idx0 = u32(_pext_u64(1, PEXT_MASK));
#else
  idx0 = 0;
  for (int i = 0; i < 16; ++i) idx0 |= u32((u64(1) >> picked[i]) & 1) << i;
#endif
  PEXT_IDX_local[0] = idx0;
 }
 u64 contribution[16] = {};
 bool contrib_found[16] = {};
 if (__builtin_popcount(PEXT_IDX_local[0]) == 1) {
  int b = __builtin_ctz(PEXT_IDX_local[0]);
  contribution[b] = 1;
  contrib_found[b] = true;
 }
 for (int k = 1; k < 65535; ++k) {
  PW_local[k] = mul(PW_local[k-1], SIGMA);
  u32 idx;
#if HAVE_PEXT
  idx = u32(_pext_u64(PW_local[k], PEXT_MASK));
#else
  idx = 0;
  for (int i = 0; i < 16; ++i) idx |= u32((PW_local[k] >> picked[i]) & 1) << i;
#endif
  PEXT_IDX_local[k] = idx;
  // single-bit idx を contribution として保存
  if (__builtin_popcount(idx) == 1) {
   int bit_pos = __builtin_ctz(idx);
   if (!contrib_found[bit_pos]) {
    contribution[bit_pos] = PW_local[k];
    contrib_found[bit_pos] = true;
   }
  }
 }
 // INV_PEXT_IDX[idx_k] = idx_{(65535 - k) % 65535}
 for (u32 idx = 0; idx < 65536; ++idx) INV_PEXT_IDX[idx] = 0;  // default
 for (int k = 0; k < 65535; ++k) {
  int inv_k = (k == 0) ? 0 : (65535 - k);
  INV_PEXT_IDX[PEXT_IDX_local[k]] = u16(PEXT_IDX_local[inv_k]);
 }
 // EMBED_BYTE 構築
 for (int p = 0; p < 2; ++p) {
  for (int b = 0; b < 256; ++b) {
   u64 v = 0;
   for (int bit = 0; bit < 8; ++bit) {
    if ((b >> bit) & 1) v ^= contribution[p * 8 + bit];
   }
   EMBED_BYTE[p][b] = v;
  }
 }
}

[[gnu::always_inline]] inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}

[[gnu::always_inline]] inline u32 extract_idx(u64 N) {
#if HAVE_PEXT
 return u32(_pext_u64(N, PEXT_MASK));
#else
 u32 r = 0;
 for (int i = 0; i < 16; ++i) r |= u32((N >> PEXT_POS[i]) & 1) << i;
 return r;
#endif
}

[[gnu::always_inline]] inline u64 embed_idx(u16 idx) {
 return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)];
}

PCLMUL_FN u64 inv_in_f16(u64 N_poly) {
 return embed_idx(INV_PEXT_IDX[extract_idx(N_poly)]);
}

PCLMUL_FN void inv_batch4(const u64 a[4], u64 out[4]) {
 u64 b1[4], b2[4], b3[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b1[k] = frob16(a[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b2[k] = frob16(b1[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b3[k] = frob16(b2[k]);
 u64 beta[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) beta[k] = mul(mul(b1[k], b2[k]), b3[k]);
 u64 N[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) N[k] = mul(a[k], beta[k]);
 u64 N_inv[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) N_inv[k] = inv_in_f16(N[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) out[k] = mul(beta[k], N_inv[k]);
}

PCLMUL_FN u64 inv_single(u64 a) {
 const u64 b1 = frob16(a);
 const u64 b2 = frob16(b1);
 const u64 b3 = frob16(b2);
 const u64 beta = mul(mul(b1, b2), b3);
 const u64 N = mul(a, beta);
 return mul(beta, inv_in_f16(N));
}

} // namespace gf2_64_pclmul_norm_pext_compact

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_pext_compact::inv_batch4;
  using gf2_64_pclmul_norm_pext_compact::inv_single;
  using gf2_64_pclmul_norm_pext_compact::init_tables;
  init_tables();
  const size_t T = as.size();
  vector<u64> ans(T);
  size_t i = 0;
  for (; i + 4 <= T; i += 4) {
   u64 b_inv[4];
   inv_batch4(&bs[i], b_inv);
   #pragma GCC unroll 4
   for (int k = 0; k < 4; ++k) ans[i+k] = mul(as[i+k], b_inv[k]);
  }
  for (; i < T; ++i) ans[i] = mul(as[i], inv_single(bs[i]));
  return ans;
 }
};
