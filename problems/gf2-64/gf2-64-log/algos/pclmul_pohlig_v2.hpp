#pragma once
// Pohlig-Hellman log_2(x) in F_{2^64}^*、Nimber.hpp 流の最適化を取り入れた v2。
//
// 主な工夫 (v1 = pclmul_pohlig.hpp 比):
//   1. F_{2^16}^* (順序 65535 = 3·5·17·257) を 1 個の subgroup として扱う。
//      4 個の素因数 BSGS の代わりに **F_{2^16} log/exp 直引き** で O(1) に。
//      → 7 BSGS → 3 BSGS + 1 log/exp。
//   2. バケットソート BSGS (Nimber 風): hash collision なし、線形探索が L1 内
//   3. CRT 段階合成の逆元は init で事前計算。
//   4. proj 用指数を hex 直書き (compile-time 定数として確実に optimize)。
//
// 期待: v1 比 ~30-40% 高速化。
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

namespace gf2_64_log_pohlig_v2 {
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

// =============================================================================
// F_{2^16}^* log/exp tables (Nimber 互換)
// =============================================================================
inline u16 PW16[65536], LN16[65536];

PCLMUL_FN void init_f16_logexp() {
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

// =============================================================================
// BSGS (open addressing hash, sqrt(p) optimal m)
// =============================================================================
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
};

inline BSGSTable bsgs_641, bsgs_65537, bsgs_6700417;

// =============================================================================
// F_{2^16}^* 内の log:
//   y = x^M ∈ F_{2^16}^*. h = g^M (g = 2 ∈ F_{2^64}^*) も F_{2^16}^* primitive。
//   y = h^{r1} ⇒ r1 = log_h(y) = log_s(y) * log_s(h)^{-1} mod 65535
//   ここで s は Nimber 流の F_{2^16}^* primitive (LN16 の基準)。
//   log_s(h) は init で計算 (Nimber のテーブルでは a = LN16[h_nim_16] という非自明な値)。
// =============================================================================
inline u32 H_LOG_INV;  // = log_s(h)^{-1} mod 65535 (init で計算)

PCLMUL_FN u32 solve_f16(u64 x_proj_poly) {
 const u64 x_nim = gf2_64_basis::poly_to_nim(x_proj_poly);
 const u16 x16 = u16(x_nim);
 if (x16 == 0) return 0;
 const u32 log_x = LN16[x16];
 return u32((u64(log_x) * H_LOG_INV) % 65535);
}

// =============================================================================
// CRT 段階合成 (Garner 形式) の事前計算逆元
// 順序: r1 (mod 65535) → r0 (mod 641) → r2 (mod 65537) → r3 (mod 6700417)
// result = r1 + 65535 * t0 + (65535*641) * t2 + (65535*641*65537) * t3
//   t0 = (r0 - r1) * 65535^{-1} mod 641
//   t2 = ((r2 - r1) - 65535 * t0) * (65535 * 641)^{-1} mod 65537
//   t3 = (((r3 - r1) - 65535 * t0) - (65535 * 641) * t2) * (65535*641*65537)^{-1} mod 6700417
// 全部 compile-time constexpr で計算
// =============================================================================
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
constexpr u64 P_F17 = 65537;  // F = Fermat (2^16+1)
constexpr u64 P_BIG = 6700417;

// 65535^{-1} mod 641
constexpr u64 INV_65535_641 = modinv_runtime(65535 % 641, 641);
// (65535 * 641)^{-1} mod 65537
constexpr u64 MOD2 = (P_F16 * P_641) % P_F17;
constexpr u64 INV_MOD2_F17 = modinv_runtime(MOD2, P_F17);
// (65535 * 641 * 65537)^{-1} mod 6700417
constexpr u64 MOD3 = (P_F16 * P_641 * P_F17) % P_BIG;
constexpr u64 INV_MOD3_BIG = modinv_runtime(MOD3, P_BIG);
// 中間 mod 値 (for CRT 中の reduce)
constexpr u64 MOD_F16 = P_F16;
constexpr u64 MOD_F16_641 = P_F16 * P_641;
constexpr u64 MOD_F16_641_F17 = P_F16 * P_641 * P_F17;

constexpr u64 N_ORDER = ~u64(0);  // 2^64 - 1
constexpr u64 G_2 = 2;

// 各 projection の指数 = (2^64-1)/p
constexpr u64 EXP_F16 = 0x0001000100010001ull;       // = N/65535 = 2^48+2^32+2^16+1
constexpr u64 EXP_641 = 0x00663d80ff99c27full;       // = N/641 ... 別途検証
constexpr u64 EXP_F17 = 0x0000ffff0000ffffull;       // = N/65537
constexpr u64 EXP_BIG = 0x00000280fffffd7full;       // = N/6700417

inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // frob4 byte table
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
 init_f16_logexp();
// F_{2^16}^* primitive h = g^M、log_s(h) を計算
 const u64 h_poly = pow_bw(G_2, EXP_F16);
 const u64 h_nim = gf2_64_basis::poly_to_nim(h_poly);
 const u16 h16 = u16(h_nim);
 const u32 h_log = LN16[h16];
 H_LOG_INV = u32(modinv_runtime(h_log, 65535));
 // BSGS tables for 641, 65537, 6700417
 const u64 g_641     = pow_bw(G_2, EXP_641);
 const u64 g_65537   = pow_bw(G_2, EXP_F17);
 const u64 g_6700417 = pow_bw(G_2, EXP_BIG);
 bsgs_641     .build(g_641,     P_641);
 bsgs_65537   .build(g_65537,   P_F17);
 bsgs_6700417 .build(g_6700417, P_BIG);
}

PCLMUL_FN u64 log_g(u64 x) {
 // 4 つの projections
 const u64 x_f16 = pow_bw(x, EXP_F16);
 const u64 x_641 = pow_bw(x, EXP_641);
 const u64 x_65537 = pow_bw(x, EXP_F17);
 const u64 x_6700417 = pow_bw(x, EXP_BIG);
 // 4 つの residue
 const u32 r1 = solve_f16(x_f16);
 const u32 r0 = bsgs_641     .solve(x_641);
 const u32 r2 = bsgs_65537   .solve(x_65537);
 const u32 r3 = bsgs_6700417 .solve(x_6700417);
 // Garner 形式 CRT
 // result = r1 + 65535 * t0 + 65535*641 * t2 + 65535*641*65537 * t3
 const u64 cur_mod_641 = u64(r1) % P_641;
 const u64 diff0 = (r0 + P_641 - cur_mod_641) % P_641;
 const u64 t0 = (diff0 * INV_65535_641) % P_641;
 //
 const u64 cur_after0_mod_F17 = (u64(r1) + MOD_F16 * t0) % P_F17;
 const u64 diff2 = (r2 + P_F17 - cur_after0_mod_F17) % P_F17;
 const u64 t2 = (diff2 * INV_MOD2_F17) % P_F17;
 //
 const u64 cur_after2_mod_BIG = (u64(r1) + (MOD_F16 * t0) % P_BIG + (MOD_F16_641 * t2) % P_BIG) % P_BIG;
 const u64 diff3 = (r3 + P_BIG - cur_after2_mod_BIG) % P_BIG;
 const u64 t3 = (__uint128_t(diff3) * INV_MOD3_BIG) % P_BIG;
 //
 return u64(r1) + MOD_F16 * t0 + MOD_F16_641 * t2 + MOD_F16_641_F17 * t3;
}

} // namespace gf2_64_log_pohlig_v2

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig_v2::log_g;
  using gf2_64_log_pohlig_v2::init_tables;
  init_tables();
  vector<u64> ans(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) ans[i] = log_g(xs[i]);
  return ans;
 }
};
