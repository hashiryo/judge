#pragma once
// pclmul_pohlig_v4.hpp の発展形:
// 6700417 BSGS を 4-way batch 化して complex memory latency を OoO で隠蔽。
//
// 単一 BSGS の giant step は (1 mul + 1 hash lookup) の連鎖で、各 ~10 cycle 直列。
// 4 query を batch して step 同期で進めれば、PCLMUL + load を並列発行できる:
//   - 4 muls per step → PCLMUL throughput 1/cycle で 4 cycle
//   - 4 lookups per step → load 2/cycle で 2 cycle
//   - 並列で ~5-6 cycle / step (vs 10 cycle 直列)
// → per-query で ~10-15K cycle 削減 (= 50% 程度短縮)
//
// implementation note:
//   - 早期終了は branch コストで OoO を阻害するので flat loop
//   - 各 query の log は m iter 全走査後にまとめて取得
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

namespace gf2_64_log_pohlig_v5 {
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

inline u64 FROB4_BYTE[8][256];

PCLMUL_FN u64 frob4(u64 a) {
 return FROB4_BYTE[0][u8(a)]       ^ FROB4_BYTE[1][u8(a >>  8)]
      ^ FROB4_BYTE[2][u8(a >> 16)] ^ FROB4_BYTE[3][u8(a >> 24)]
      ^ FROB4_BYTE[4][u8(a >> 32)] ^ FROB4_BYTE[5][u8(a >> 40)]
      ^ FROB4_BYTE[6][u8(a >> 48)] ^ FROB4_BYTE[7][u8(a >> 56)];
}

PCLMUL_FN u64 pow_bw(u64 a, u64 e) {
 if (e == 0) return 1;
 u64 T[16];
 T[0] = 1; T[1] = a;
 #pragma GCC unroll 14
 for (int i = 2; i < 16; ++i) T[i] = mul(T[i-1], a);
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

constexpr u64 SIGMA = 0xa1573a4da2bc3a32ull;
inline u64 PEXT_MASK = 0;
inline u16 LN_SIGMA_PEXT[65536];
#if !HAVE_PEXT
inline int PEXT_POS[16];
#endif

PCLMUL_FN void build_sigma_pext_table() {
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
 u64 cur = 1;
 for (u32 k = 0; k < 65535; ++k) {
  u32 idx;
#if HAVE_PEXT
  idx = u32(_pext_u64(cur, PEXT_MASK));
#else
  idx = 0;
  for (int i = 0; i < 16; ++i) idx |= u32((cur >> picked[i]) & 1) << i;
#endif
  LN_SIGMA_PEXT[idx] = u16(k);
  cur = mul(cur, SIGMA);
 }
 LN_SIGMA_PEXT[0] = 0;
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

struct DirectLogTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u64 mask;

 PCLMUL_RUN void build(u64 base, u32 p) {
  u64 cap = 8;
  while (cap < 2u * p) cap *= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask = cap - 1;
  u64 cur = 1;
  for (u32 k = 0; k < p; ++k) {
   u64 h = (cur * 0x9E3779B97F4A7C15ull) & mask;
   while (keys[h] != ~u64(0)) h = (h + 1) & mask;
   keys[h] = cur;
   values[h] = k;
   cur = mul(cur, base);
  }
 }

 PCLMUL_RUN u32 lookup(u64 target) const {
  u64 h = (target * 0x9E3779B97F4A7C15ull) & mask;
  while (keys[h] != ~u64(0)) {
   if (keys[h] == target) return values[h];
   h = (h + 1) & mask;
  }
  return ~u32(0);
 }
};

inline DirectLogTable direct_641, direct_65537;

struct BSGSTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u64 mask;
 u64 m;
 u64 q;
 u64 inv_base_m;

 PCLMUL_RUN void build(u64 base, u64 q_) {
  q = q_;
  m = 1;
  while (m * m < q) ++m;
  u64 cap = 8;
  while (cap < 4 * m) cap *= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask = cap - 1;
  u64 cur = 1;
  for (u64 j = 0; j < m; ++j) {
   u64 h = (cur * 0x9E3779B97F4A7C15ull) & mask;
   while (keys[h] != ~u64(0)) h = (h + 1) & mask;
   keys[h] = cur;
   values[h] = u32(j);
   cur = mul(cur, base);
  }
  inv_base_m = pow_bw(base, q - m);
 }

 PCLMUL_RUN u32 solve(u64 target) const {
  u64 t = target;
  for (u64 i = 0; i < m; ++i) {
   u64 h = (t * 0x9E3779B97F4A7C15ull) & mask;
   while (keys[h] != ~u64(0)) {
    if (keys[h] == t) {
     u64 res = i * m + values[h];
     if (res < q) return u32(res);
    }
    h = (h + 1) & mask;
   }
   t = mul(t, inv_base_m);
  }
  return u32(q);
 }

 // 4-way batch BSGS: 4 query を同期した giant step で並列処理
 // OoO で 4 つの (mul + lookup) を並行発行、memory latency 隠蔽
 PCLMUL_RUN void solve_batch4(const u64 targets[4], u32 results[4]) const {
  u64 t[4]   = {targets[0], targets[1], targets[2], targets[3]};
  u32 res[4] = {u32(q), u32(q), u32(q), u32(q)};
  for (u64 i = 0; i < m; ++i) {
   // 4 lookups 並列 (PCLMUL とは異なる port、並行発行可)
   #pragma GCC unroll 4
   for (int k = 0; k < 4; ++k) {
    u64 h = (t[k] * 0x9E3779B97F4A7C15ull) & mask;
    while (keys[h] != ~u64(0)) {
     if (keys[h] == t[k] && res[k] == u32(q)) {
      u64 r = i * m + values[h];
      if (r < q) res[k] = u32(r);
     }
     h = (h + 1) & mask;
    }
   }
   // 4 muls 並列 (互いに依存なし、PCLMUL 1/cycle throughput で 4 cycle)
   #pragma GCC unroll 4
   for (int k = 0; k < 4; ++k) t[k] = mul(t[k], inv_base_m);
  }
  for (int k = 0; k < 4; ++k) results[k] = res[k];
 }
};

inline BSGSTable bsgs_6700417;
inline u32 H_LOG_INV;

PCLMUL_FN u32 solve_f16(u64 x_proj_poly) {
 const u32 idx = extract_idx(x_proj_poly);
 if (idx == 0) return 0;
 const u32 log_x = LN_SIGMA_PEXT[idx];
 return u32((u64(log_x) * H_LOG_INV) % 65535);
}

constexpr u64 modinv_runtime(u64 a, u64 m) {
 long long r0 = (long long) m, r1 = (long long) a;
 long long s0 = 0, s1 = 1;
 while (r1 != 0) {
  long long q = r0 / r1;
  long long r2 = r0 - q * r1; r0 = r1; r1 = r2;
  long long s2 = s0 - q * s1; s0 = s1; s1 = s2;
 }
 if (s0 < 0) s0 += (long long) m;
 return (u64) s0;
}

constexpr u64 P_F16 = 65535;
constexpr u64 P_641 = 641;
constexpr u64 P_F17 = 65537;
constexpr u64 P_BIG = 6700417;

constexpr u64 INV_65535_641 = modinv_runtime(65535 % 641, 641);
constexpr u64 MOD2 = (P_F16 * P_641) % P_F17;
constexpr u64 INV_MOD2_F17 = modinv_runtime(MOD2, P_F17);
constexpr u64 MOD3 = (P_F16 * P_641 * P_F17) % P_BIG;
constexpr u64 INV_MOD3_BIG = modinv_runtime(MOD3, P_BIG);
constexpr u64 MOD_F16 = P_F16;
constexpr u64 MOD_F16_641 = P_F16 * P_641;
constexpr u64 MOD_F16_641_F17 = P_F16 * P_641 * P_F17;

constexpr u64 G_2 = 2;
constexpr u64 EXP_F16 = 0x0001000100010001ull;
constexpr u64 EXP_641 = 0x00663d80ff99c27full;
constexpr u64 EXP_F17 = 0x0000ffff0000ffffull;
constexpr u64 EXP_BIG = 0x00000280fffffd7full;

inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
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
 build_sigma_pext_table();
 const u64 h_poly = pow_bw(G_2, EXP_F16);
 const u32 h_idx = extract_idx(h_poly);
 const u32 h_log = LN_SIGMA_PEXT[h_idx];
 H_LOG_INV = u32(modinv_runtime(h_log, 65535));
 const u64 g_641     = pow_bw(G_2, EXP_641);
 const u64 g_65537   = pow_bw(G_2, EXP_F17);
 const u64 g_6700417 = pow_bw(G_2, EXP_BIG);
 direct_641   .build(g_641,     u32(P_641));
 direct_65537 .build(g_65537,   u32(P_F17));
 bsgs_6700417 .build(g_6700417, P_BIG);
}

PCLMUL_FN u64 combine_crt(u32 r1, u32 r0, u32 r2, u32 r3) {
 const u64 cur_mod_641 = u64(r1) % P_641;
 const u64 diff0 = (r0 + P_641 - cur_mod_641) % P_641;
 const u64 t0 = (diff0 * INV_65535_641) % P_641;
 const u64 cur_after0_mod_F17 = (u64(r1) + MOD_F16 * t0) % P_F17;
 const u64 diff2 = (r2 + P_F17 - cur_after0_mod_F17) % P_F17;
 const u64 t2 = (diff2 * INV_MOD2_F17) % P_F17;
 const u64 cur_after2_mod_BIG = (u64(r1) + (MOD_F16 * t0) % P_BIG + (MOD_F16_641 * t2) % P_BIG) % P_BIG;
 const u64 diff3 = (r3 + P_BIG - cur_after2_mod_BIG) % P_BIG;
 const u64 t3 = (__uint128_t(diff3) * INV_MOD3_BIG) % P_BIG;
 return u64(r1) + MOD_F16 * t0 + MOD_F16_641 * t2 + MOD_F16_641_F17 * t3;
}

// 4 query をまとめて log_g、6700417 BSGS を batch 化
PCLMUL_FN void log_g_batch4(const u64 xs[4], u64 results[4]) {
 u64 x_f16[4], x_641[4], x_65537[4], x_6700417[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) {
  x_f16[k]     = pow_bw(xs[k], EXP_F16);
  x_641[k]     = pow_bw(xs[k], EXP_641);
  x_65537[k]   = pow_bw(xs[k], EXP_F17);
  x_6700417[k] = pow_bw(xs[k], EXP_BIG);
 }
 u32 r1[4], r0[4], r2[4], r3[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) {
  r1[k] = solve_f16(x_f16[k]);
  r0[k] = direct_641   .lookup(x_641[k]);
  r2[k] = direct_65537 .lookup(x_65537[k]);
 }
 // ここが本命: 6700417 BSGS を batch4 で
 bsgs_6700417.solve_batch4(x_6700417, r3);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) results[k] = combine_crt(r1[k], r0[k], r2[k], r3[k]);
}

// Single query (T が 4 の倍数でない時の余り処理用)
PCLMUL_FN u64 log_g_single(u64 x) {
 const u64 x_f16 = pow_bw(x, EXP_F16);
 const u64 x_641 = pow_bw(x, EXP_641);
 const u64 x_65537 = pow_bw(x, EXP_F17);
 const u64 x_6700417 = pow_bw(x, EXP_BIG);
 const u32 r1 = solve_f16(x_f16);
 const u32 r0 = direct_641   .lookup(x_641);
 const u32 r2 = direct_65537 .lookup(x_65537);
 const u32 r3 = bsgs_6700417 .solve(x_6700417);
 return combine_crt(r1, r0, r2, r3);
}

} // namespace gf2_64_log_pohlig_v5

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig_v5::log_g_batch4;
  using gf2_64_log_pohlig_v5::log_g_single;
  using gf2_64_log_pohlig_v5::init_tables;
  init_tables();
  const size_t T = xs.size();
  vector<u64> ans(T);
  size_t i = 0;
  for (; i + 4 <= T; i += 4) {
   u64 ys[4] = {xs[i], xs[i+1], xs[i+2], xs[i+3]};
   u64 rs[4];
   log_g_batch4(ys, rs);
   ans[i] = rs[0]; ans[i+1] = rs[1]; ans[i+2] = rs[2]; ans[i+3] = rs[3];
  }
  for (; i < T; ++i) ans[i] = log_g_single(xs[i]);
  return ans;
 }
};
