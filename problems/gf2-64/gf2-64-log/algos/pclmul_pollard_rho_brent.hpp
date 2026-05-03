#pragma once
// =============================================================================
// 6700417 subgroup の log を Pollard rho + Brent's cycle detection で計算する版。
// pclmul_pohlig_v5.hpp の batch4 BSGS と置換可能な代替実装。
//
// 結論: BSGS より遅い (実測 random_00 で 1900ms vs v5 390ms ≈ 5× 遅い)。
//   pohlig_v5 を超えない。教育/比較用に保持。
//
// 速度劣化の主な要因:
//   1. per-step の mod operation: a, b ∈ Z/6700417 を毎 step 更新するため
//      ~5 cycle × 2 mod = 10 cycle、BSGS の hash lookup ~5 cycle より重い
//   2. Brent's expected iter ~3245 だが分散大、最悪 BSGS の m=2589 を超える
//   3. memory が定数 (=hash table 不要) という Pollard rho の利点が、
//      BSGS table が L2 (~200KB) で済む我々の状況では無効
//
// Pollard rho が勝つのは:
//   - p が極端に大きく BSGS table が cache に乗らない場合 (例: p ≥ 10^12)
//   - memory-constrained 環境
//   - Pollard kangaroo で SIMD 多並列化する研究シナリオ
//
// 我々の p = 6700417 では BSGS が一方的に強い、というのが実証された。
// =============================================================================
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

namespace gf2_64_log_pollard_rho_brent {
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

// =============================================================================
// Pollard rho with Brent's cycle detection for 6700417 subgroup
// =============================================================================
constexpr u32 P_BIG = 6700417;
inline u64 g_6700417;  // base = g^{N/6700417}

// mod p の逆元 (Fermat: a^{p-2})
[[gnu::always_inline]] inline u32 mod_inv_p(u32 a, u32 p) {
 u32 r = 1, b = a, e = p - 2;
 while (e) { if (e & 1) r = u64(r) * b % p; b = u64(b) * b % p; e >>= 1; }
 return r;
}

// Pollard rho walk function: state = (α, a, b) where α = base^a · target^b (in subgroup)
// 4 partitions based on α & 0x3:
//   0, 1: multiply by base   (a += 1)
//   2:    square              (a, b doubled)
//   3:    multiply by target  (b += 1)
struct WalkState { u64 alpha; u32 a; u32 b; };

[[gnu::always_inline]] inline WalkState rho_step(WalkState s, u64 base, u64 target) {
 const u32 sel = u32(s.alpha & 0x3);
 if (sel < 2) {
  s.alpha = mul(s.alpha, base);
  s.a = s.a + 1; if (s.a >= P_BIG) s.a -= P_BIG;
 } else if (sel == 2) {
  s.alpha = sq(s.alpha);
  s.a = u64(s.a) * 2 % P_BIG;
  s.b = u64(s.b) * 2 % P_BIG;
 } else {
  s.alpha = mul(s.alpha, target);
  s.b = s.b + 1; if (s.b >= P_BIG) s.b -= P_BIG;
 }
 return s;
}

// Pollard rho + Brent's cycle detection
PCLMUL_FN u32 pollard_rho_log(u64 base, u64 target) {
 if (target == 1) return 0;
 // 初期状態を (target, 0, 1) にすることで b > 0 を保証 (collision で b 差が出やすい)
 WalkState slow = {target, 0, 1};
 WalkState fast = rho_step(slow, base, target);
 u32 power = 1, lam = 1;
 // 安全上限: Brent's expected ~3245、保険で 10× = 32K iter で諦める
 constexpr u32 MAX_ITER = 65536;
 for (u32 it = 0; it < MAX_ITER && slow.alpha != fast.alpha; ++it) {
  if (power == lam) {
   slow = fast;
   power *= 2;
   lam = 0;
  }
  fast = rho_step(fast, base, target);
  ++lam;
 }
 if (slow.alpha != fast.alpha) return u32(-1);  // failure
 // 衝突: base^slow.a · target^slow.b = base^fast.a · target^fast.b
 // ⇒ base^{slow.a - fast.a} = target^{fast.b - slow.b}
 // ⇒ x · (fast.b - slow.b) ≡ slow.a - fast.a mod P_BIG
 // ⇒ x = (slow.a - fast.a) · (fast.b - slow.b)^{-1} mod P_BIG
 u32 db = fast.b >= slow.b ? fast.b - slow.b : P_BIG - (slow.b - fast.b);
 u32 da = slow.a >= fast.a ? slow.a - fast.a : P_BIG - (fast.a - slow.a);
 if (db == 0) return u32(-1);  // rare, restart needed (we just fail)
 u32 inv_db = mod_inv_p(db, P_BIG);
 return u64(da) * inv_db % P_BIG;
}

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
 g_6700417 = pow_bw(G_2, EXP_BIG);
 direct_641   .build(pow_bw(G_2, EXP_641), u32(P_641));
 direct_65537 .build(pow_bw(G_2, EXP_F17), u32(P_F17));
}

PCLMUL_FN u64 log_g(u64 x) {
 const u64 x_f16 = pow_bw(x, EXP_F16);
 const u64 x_641 = pow_bw(x, EXP_641);
 const u64 x_65537 = pow_bw(x, EXP_F17);
 const u64 x_6700417 = pow_bw(x, EXP_BIG);
 const u32 r1 = solve_f16(x_f16);
 const u32 r0 = direct_641   .lookup(x_641);
 const u32 r2 = direct_65537 .lookup(x_65537);
 const u32 r3 = pollard_rho_log(g_6700417, x_6700417);
 // Garner CRT
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

} // namespace gf2_64_log_pollard_rho_brent

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pollard_rho_brent::log_g;
  using gf2_64_log_pollard_rho_brent::init_tables;
  init_tables();
  vector<u64> ans(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) ans[i] = log_g(xs[i]);
  return ans;
 }
};
