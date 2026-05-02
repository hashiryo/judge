#pragma once
// Norm-based + PEXT-based subfield inv (Nimber 依存ゼロ)。
//
// 要点:
//   F_{2^16} ⊂ F_{2^64} は σ^0, σ^1, ..., σ^{15} の張る 16 次元 GF(2) 部分空間。
//   この部分空間に「適切な 16 bit 位置」での射影が全単射になる、というプロパティ
//   が常に存在する (= 16 個の線型独立な座標を選べる)。
//   それらの bit を MASK で集めて PEXT 命令で 1 サイクルで 16-bit 表現を抽出。
//
//   そして INV_PEXT[16-bit index] = 64-bit poly の inverse  を init で構築すれば、
//   subfield inv は **PEXT 1 + lookup 1 = 2 命令** で完了。
//
// メモリ: INV_PEXT = 65536 × 8 = 512 KiB (L2 hit)
// 命令数: pclmul_norm_inv16_compact (8+1+2 = 11 lookups) → 1 lookup + 1 PEXT
// ARM: PEXT 無いので _shared 経由 fallback (= 普通の subfield inv)
//
// 注: σ は gen_tower_custom_sparse から取った 0xa1573a4da2bc3a32
//     (r(s) = s^16 + s^12 + s^3 + s + 1 で r(σ) = 0)。
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

namespace gf2_64_pclmul_norm_pext {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;
using gf2_64_pclmul::sq;
using u16 = unsigned short;
using u32 = unsigned;

// F_{2^16} subfield 生成元 (poly basis、r(σ) = 0 で r = s^16 + s^12 + s^3 + s + 1)
constexpr u64 SIGMA = 0xa1573a4da2bc3a32ull;

inline u64 FROB16_BYTE[8][256];
inline u64 INV_PEXT[65536];   // 16-bit PEXT idx → 64-bit poly inverse
inline u64 PEXT_MASK = 0;     // 16 個の線型独立 bit 位置の mask

#if !HAVE_PEXT
// ARM fallback: 16 個の bit 位置から手動抽出するための shift/and 列を保持
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

 // σ^0..σ^{15} を計算
 u64 sigma_pow[16];
 sigma_pow[0] = 1;
 for (int i = 1; i < 16; ++i) sigma_pow[i] = mul(sigma_pow[i-1], SIGMA);

 // 64×16 行列 M_ij = bit i of sigma_pow[j] から、線型独立な 16 個の row を選ぶ
 // (Gauss 消去でピボット行を見つける)
 u16 row_vec[64];  // 各 row の 16-bit 表現
 for (int r = 0; r < 64; ++r) {
  u16 v = 0;
  for (int c = 0; c < 16; ++c) {
   if ((sigma_pow[c] >> r) & 1) v |= u16(1) << c;
  }
  row_vec[r] = v;
 }
 // 線型独立な 16 行を選ぶ (各 basis[k] は "k 番目の bit が leading" を保持)
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
 // PEXT_MASK = picked bits の OR
 PEXT_MASK = 0;
 for (int i = 0; i < 16; ++i) PEXT_MASK |= (u64(1) << picked[i]);
#if !HAVE_PEXT
 for (int i = 0; i < 16; ++i) PEXT_POS[i] = picked[i];
#endif

 // 「PEXT idx → 元の 16 σ-coeff」の変換行列を準備
 // PEXT(N, MASK) は picked 順に bit を集めるので、 i 番目の出力 bit = N の bit picked[i]
 // 一方「σ-coeff」は σ^i 軸上の係数。σ^i の picked[k] bit = (sigma_pow[i] >> picked[k]) & 1
 // 16×16 行列 M: row k = picked[k] bit、col i = σ^i 軸 → M[k][i] = (sigma_pow[i] >> picked[k]) & 1
 // PEXT 出力 = M · (σ-coeff vec)
 // 逆: σ-coeff = M^{-1} · PEXT 出力
 u32 M[16];      // M[k] = row k as 16-bit
 u32 Minv[16];   // M^{-1}
 for (int k = 0; k < 16; ++k) {
  u32 v = 0;
  for (int i = 0; i < 16; ++i) {
   if ((sigma_pow[i] >> picked[k]) & 1) v |= u32(1) << i;
  }
  M[k] = v;
  Minv[k] = u32(1) << k;
 }
 // M を単位行列に変換しつつ Minv を構築 (拡大行列で Gauss-Jordan)
 for (int col = 0; col < 16; ++col) {
  // pivot を col 行に持ってくる
  int pv = -1;
  for (int r = col; r < 16; ++r) if ((M[r] >> col) & 1) { pv = r; break; }
  // ありえないが念の為
  if (pv == -1) continue;
  if (pv != col) { std::swap(M[col], M[pv]); std::swap(Minv[col], Minv[pv]); }
  // 他の行から col 列を消去
  for (int r = 0; r < 16; ++r) {
   if (r != col && ((M[r] >> col) & 1)) { M[r] ^= M[col]; Minv[r] ^= Minv[col]; }
  }
 }
 // 今 Minv は σ-coeff = Minv · PEXT_idx の行列

 // INV_PEXT[idx] を構築:
 //   1) σ-coeff c = Minv · idx
 //   2) N (poly 64-bit) = ⊕ c_i σ^i
 //   3) inv_N = N^{2^16 - 2}  (F_{2^16} 内の inv、Itoh-Tsujii で計算)
 //   4) INV_PEXT[idx] = inv_N
 // INV_PEXT を σ^k 列挙で構築 (~30 μs):
 //   σ^k for k = 0..65534 を 65534 PCLMUL muls で chain 計算
 //   各々の PEXT idx と σ^{-k} を対応付け
 // 0 は subfield に住まないので INV_PEXT[pext(0)] = 0 とだけ設定 (= INV_PEXT[0])
 for (u32 idx = 0; idx < 65536; ++idx) INV_PEXT[idx] = 0;  // 初期化
 // σ^k を chain で生成、同時に inv pair を埋める
 // σ^{-k} = σ^{65535 - k} (k=0 のとき σ^0 = 1 で自分自身 inv)
 // 一旦 PW[k] = σ^k を保持してから inv ペアを埋める (= σ^k と σ^{(65535-k)%65535})
 // PW_local は init 後に解放したいので vector
 std::vector<u64> PW_local(65535);
 PW_local[0] = 1;
 for (int k = 1; k < 65535; ++k) {
  PW_local[k] = mul(PW_local[k-1], SIGMA);
 }
 for (int k = 0; k < 65535; ++k) {
  u32 pext_k;
#if HAVE_PEXT
  pext_k = u32(_pext_u64(PW_local[k], PEXT_MASK));
#else
  pext_k = 0;
  for (int i = 0; i < 16; ++i) {
   pext_k |= u32((PW_local[k] >> picked[i]) & 1) << i;
  }
#endif
  int inv_k = (k == 0) ? 0 : (65535 - k);
  INV_PEXT[pext_k] = PW_local[inv_k];
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
 for (int i = 0; i < 16; ++i) {
  r |= u32((N >> PEXT_POS[i]) & 1) << i;
 }
 return r;
#endif
}

PCLMUL_FN u64 inv_in_f16(u64 N_poly) {
 return INV_PEXT[extract_idx(N_poly)];
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

} // namespace gf2_64_pclmul_norm_pext

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_pext::inv_batch4;
  using gf2_64_pclmul_norm_pext::inv_single;
  using gf2_64_pclmul_norm_pext::init_tables;
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
