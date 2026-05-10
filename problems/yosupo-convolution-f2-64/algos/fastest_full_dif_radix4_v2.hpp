#pragma once
#include "_common.hpp"
// fastest_full_dif_radix4.hpp の **v2** 版。
// 違い: q=1 (level pair (2, 1)) の scalar fallback も VPCLMULQDQ 化。
//
// 元の scalar 内訳 (4 個独立 mul):
//   Step 1: mul(Q2, v), mul(Q3, v)              ← 同 scalar v
//   Step 2: mul(A1, u_lo), mul(A3, u_hi)        ← 異なる scalar
//
// Step 1: pack (Q2, Q3) + broadcast(v) → 1 VPCLMULQDQ
// Step 2: pack (A1, A3) + pack(u_lo, u_hi)     → 1 VPCLMULQDQ
//   (pack_u_pair で lane 0 = u_lo, lane 2 = u_hi の vector を 1 度作って block 内で再利用)
//
// 期待効果: q=1 level (= 全 mul の ~1/d ≈ 5%) で 4 PCLMUL → 2 VPCLMULQDQ。
// 全体 ~2.5% 改善見込み。
//
// j=0 block (v=0, u_lo=0) は元の scalar fallback 維持 (mul skip 効果保つため)。
//
// 元 v1 (scalar 全部) のコメント:
//
// アルゴリズム:
//   - Phase A (true DIF, no shuffle): monomial → LCH 自然順
//     再帰的 P = Q_0 + s_{k-1}(x) · Q_1 分解 (Cantor subspace poly)。
//     Cantor で s_k(x) = Σ_{l ⊆ k bits} x^{2^l}。
//   - Phase B (DIF top-down butterflies, single ska per block):
//     LCH natural → eval natural。 ska は level 非依存な master 表 M[j] から。
//
// DIT 全部入り (fastest_full_v2) との対応:
//
// 1. Cantor 対称性 (DIT 用 2 mul → 1 mul):
//    DIF butterfly は元から (lo + ska·hi, lo + (ska+1)·hi) = (t, t ^ hi) で 1 mul。
//    DIT で言う Cantor 対称性を構造的に吸収済み。 移植不要。
//
// 2. halftable (twiddle メモリ半分):
//    DIF master M[j] = Σ_{L: j_L=1} β_{L+1} は元から size 2^(d-1)、 level 非依存
//    で 1 本だけ。 DIT 全部入りの per-level halftable よりさらにコンパクト。
//
// 3. ska=0 ブロック完全 skip + mul2 ペアリング:
//    block j=0 は ska=0 → ブロック全体が butterfly_0 (PCLMUL 無し) で済む。
//    j>0 のみ btf_fft_block_dif で SIMD 処理。
//
// 4. butterfly pair の SIMD pipeline (btf_fft_block_dif):
//    block 内の (i, i+1) ペアを __m256i で。 ska は scalar、 一度 broadcast して
//    block 全体で再利用 (DIT は g0[i] が i 毎に変わるのでロードが必要だったが
//    DIF はロード 1 回のみ)。 それだけ register pressure が低く pipeline がきれい。
//
// 5. Phase A 最終段 fusion:
//    DIF Phase A の総 op 数は DIT m-loop の ~1/5 (Σpopcount(k-1) vs Σ(k-1))。
//    levels 2, 3 (len=4, 8) は inner ops が単純 (popcount(1)=1, popcount(2)=1)、
//    特殊化して unroll する。
//
// 6. init() 簡略化: 既達。 master 1 本だけ。
//
// 追加最適化:
// 7. submask 表 (SUBMASK_TBL) を constexpr で .rodata に焼く。 bc_to_lch /
//    bc_to_mono の per-level submask 列挙を runtime build から table lookup に。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_full_dif_radix4_v2 {

using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// constexpr GF(2^64) mul / sq
// =============================================================================
constexpr void clmul128_ce(u64 a, u64 b, u64& lo_out, u64& hi_out) {
 u64 Tlo[16]= {0}, Thi[16]= {0};
 Tlo[1]= a;
 for(int v= 2; v < 16; ++v) {
  u64 plo= Tlo[v >> 1], phi= Thi[v >> 1];
  u64 nlo= plo << 1;
  u64 nhi= (phi << 1) | (plo >> 63);
  if(v & 1) nlo^= a;
  Tlo[v]= nlo;
  Thi[v]= nhi;
 }
 u64 lo= 0, hi= 0;
 for(int s= 60; s >= 0; s-= 4) {
  u64 nhi= (hi << 4) | (lo >> 60);
  u64 nlo= lo << 4;
  u32 nib= u32((b >> s) & 0xF);
  lo= nlo ^ Tlo[nib];
  hi= nhi ^ Thi[nib];
 }
 lo_out= lo;
 hi_out= hi;
}
constexpr u64 mul_ce(u64 a, u64 b) {
 u64 lo, hi;
 clmul128_ce(a, b, lo, hi);
 u64 f1l, f1h, f2l, f2h;
 clmul128_ce(hi, 0x1B, f1l, f1h);
 clmul128_ce(f1h, 0x1B, f2l, f2h);
 return lo ^ f1l ^ f2l;
}
constexpr u64 sq_ce(u64 a) { return mul_ce(a, a); }

// =============================================================================
// Cantor chain。 β_l = CHAIN[62-l]、 β_l^2 + β_l = β_{l-1}、 β_0 = 1
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// Submask 表 (proper submasks of k for k=0..MAX_SUB_IDX)
//   bc_to_lch / bc_to_mono の inner で使う x^{2^l} の oset (= 1<<l) を per-k で
//   precompute。 sub_idx = level - 1 ≤ d - 1 ≤ 62、 popcount ≤ 5 程度なので
//   per-k に 32 entries 確保すれば十分。 .rodata に焼いて runtime 構築排除。
// =============================================================================
// n は int サイズ ≤ 2^31 - 1、 d = msb(n) ≤ 30、 sub_idx = d - 1 ≤ 29 で十分。
// k ≤ 30 で popcount(k) max = 4 (popcount(15)=popcount(23)=popcount(30)=4 等)、
// proper submask 数 = 2^4 - 1 = 15 → MAX_SUBS_PER_K = 16 で足りる。
constexpr int MAX_SUB_IDX= 30;
constexpr int MAX_SUBS_PER_K= 16;
struct SubmaskTable {
 array<array<int, MAX_SUBS_PER_K>, MAX_SUB_IDX + 1> subs;
 array<int, MAX_SUB_IDX + 1> sub_n;
};
constexpr SubmaskTable build_submask_tbl() {
 SubmaskTable t{};
 for(int k= 1; k <= MAX_SUB_IDX; ++k) {
  int n= 0;
  int s= k;
  while(true) {
   s= (s - 1) & k;
   t.subs[k][n++]= 1 << s;
   if(s == 0) break;
  }
  t.sub_n[k]= n;
 }
 return t;
}
constexpr auto SUBMASK_TBL= build_submask_tbl();

// =============================================================================
// DIF master twiddle (single contiguous, level 非依存)
//   M[j] = Σ_{L: j_L=1} β_{L+1}   (j ∈ [0, 2^(d-1)))
// =============================================================================
struct nim_fft_data {
 std::vector<u64> master;
 int max_d= -1;
 void init(int d) {
  if(d <= max_d) return;
  max_d= d;
  int sz= (d <= 1) ? 1 : (1 << (d - 1));
  master.assign(sz, 0);
  for(int j= 1; j < sz; ++j) {
   int L= __builtin_ctz(j);
   master[j]= master[j ^ (1 << L)] ^ CHAIN[61 - L];
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// SIMD primitives (fastest_full_v2 から借用、 ska broadcast 化)
// =============================================================================
VPCLMUL inline __m256i expand_pair(const u64* p) {
 __m128i v= _mm_loadu_si128((const __m128i*)p);
 return _mm256_permute4x64_epi64(_mm256_castsi128_si256(v), _MM_SHUFFLE(3, 1, 2, 0));
}
VPCLMUL inline __m128i pack_pair(__m256i v) {
 return _mm256_castsi256_si128(_mm256_permute4x64_epi64(v, _MM_SHUFFLE(3, 2, 2, 0)));
}
inline const __m256i RED_TABLE_DIF_R4V2= _mm256_setr_epi8(
 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);

VPCLMUL inline __m256i clmul_reduce_pair(__m256i a_vec, __m256i b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE_DIF_R4V2, indices);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
}

// DIF FFT butterfly (block 内 ska 一定): pair (i, i+1) を SIMD で。
//   r = mul(hi, ska); lo' = lo ^ r; hi' = lo' ^ hi
// ska_vec は呼び出し側で 1 度だけ broadcast されたもの。
VPCLMUL inline void btf_fft_pair_dif(u64* f0_ptr, u64* f1_ptr, __m256i ska_vec) {
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i p_vec= clmul_reduce_pair(hi_vec, ska_vec);
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i t_vec= _mm256_xor_si256(lo_vec, p_vec);
 __m256i f1_new= _mm256_xor_si256(t_vec, hi_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(t_vec));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(f1_new));
}

// DIF IFFT butterfly: 順序逆。
//   hi_new = lo ^ hi (= s); lo_new = lo ^ s*ska
VPCLMUL inline void btf_ifft_pair_dif(u64* f0_ptr, u64* f1_ptr, __m256i ska_vec) {
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i s_vec= _mm256_xor_si256(lo_vec, hi_vec);
 __m256i p_vec= clmul_reduce_pair(s_vec, ska_vec);
 __m256i lo_new= _mm256_xor_si256(lo_vec, p_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(lo_new));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(s_vec));
}

// ska scalar → __m256i 形式 (lane 0, 2 のみ valid、 pair clmul の入力レイアウト)
// 構造: (ska, 0, ska, 0)
VPCLMUL inline __m256i broadcast_ska(u64 ska) {
 __m128i lo= _mm_set_epi64x(0, ska);
 __m256i v= _mm256_castsi128_si256(lo);
 return _mm256_inserti128_si256(v, lo, 1);
}

// 2 つの scalar を lane 0, 2 に並べた __m256i を作る (異なる scalar の SIMD pair 用)
// 構造: (a, 0, b, 0)
VPCLMUL inline __m256i pack_u_pair(u64 a, u64 b) {
 return _mm256_set_epi64x(0, b, 0, a);
}

// =============================================================================
// Phase A (true DIF): bc_to_lch / bc_to_mono
//   level 2, 3 は popcount(level-1) = 1 で単純 (1 submask {0})、 inner ops 1/i。
//   level k ≥ 4 は popcount(k-1) ≥ 1 (一般に複数 submask)。
// =============================================================================
inline void bc_to_lch(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= d; level >= 2; --level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  const int sub_n= SUBMASK_TBL.sub_n[sub_idx];
  const int* subs= SUBMASK_TBL.subs[sub_idx].data();
  if(sub_n == 1) {
   // 単一 proper submask {0}、 inner op は p[base+i+1] ^= p[base+half+i] のみ
   for(int base= 0; base < n; base+= len) {
    u64* p= poly + base;
    for(int i= half - 1; i >= 0; --i) p[i + 1]^= p[half + i];
   }
  } else {
   for(int base= 0; base < n; base+= len) {
    u64* p= poly + base;
    for(int i= half - 1; i >= 0; --i) {
     u64 q= p[half + i];
     if(q == 0) continue;
     for(int idx= 0; idx < sub_n; ++idx) p[i + subs[idx]]^= q;
    }
   }
  }
 }
}

inline void bc_to_mono(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= 2; level <= d; ++level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  const int sub_n= SUBMASK_TBL.sub_n[sub_idx];
  const int* subs= SUBMASK_TBL.subs[sub_idx].data();
  if(sub_n == 1) {
   for(int base= 0; base < n; base+= len) {
    u64* p= poly + base;
    for(int i= 0; i < half; ++i) p[i + 1]^= p[half + i];
   }
  } else {
   for(int base= 0; base < n; base+= len) {
    u64* p= poly + base;
    for(int i= 0; i < half; ++i) {
     u64 q= p[half + i];
     if(q == 0) continue;
     // forward の逆順 (subs を逆順で適用)
     for(int idx= sub_n - 1; idx >= 0; --idx) p[i + subs[idx]]^= q;
    }
   }
  }
 }
}

// =============================================================================
// Phase B: DIF top-down butterflies、 j=0 完全 skip + i pair SIMD
// =============================================================================
// radix-4 forward butterfly group (q ≥ 2、 i pair SIMD):
//   Q_j を __m256i pair で load → step1 (v) → step2 lower (u_lo) / upper (u_hi)
//   ska=0 早期分岐は呼び出し側で処理 (関数 4 つに分岐、 ここでは general case のみ)
VPCLMUL inline void btf_fft_radix4_pair(u64* base, int q, int q2, int q3, int i, __m256i v_vec, __m256i u_lo_vec, __m256i u_hi_vec) {
 __m256i Q0= expand_pair(base + i);
 __m256i Q1= expand_pair(base + q + i);
 __m256i Q2= expand_pair(base + q2 + i);
 __m256i Q3= expand_pair(base + q3 + i);
 // Step 1: A_0 = Q_0 ^ v Q_2; A_1 = Q_1 ^ v Q_3
 __m256i A0= _mm256_xor_si256(Q0, clmul_reduce_pair(Q2, v_vec));
 __m256i A1= _mm256_xor_si256(Q1, clmul_reduce_pair(Q3, v_vec));
 __m256i A2= _mm256_xor_si256(A0, Q2);
 __m256i A3= _mm256_xor_si256(A1, Q3);
 // Step 2 lower: R_C00 = A_0 ^ u_lo A_1; R_C01 = R_C00 ^ A_1
 __m256i R0= _mm256_xor_si256(A0, clmul_reduce_pair(A1, u_lo_vec));
 __m256i R1= _mm256_xor_si256(R0, A1);
 // Step 2 upper: R_C10 = A_2 ^ u_hi A_3; R_C11 = R_C10 ^ A_3
 __m256i R2= _mm256_xor_si256(A2, clmul_reduce_pair(A3, u_hi_vec));
 __m256i R3= _mm256_xor_si256(R2, A3);
 _mm_storeu_si128((__m128i*)(base + i), pack_pair(R0));
 _mm_storeu_si128((__m128i*)(base + q + i), pack_pair(R1));
 _mm_storeu_si128((__m128i*)(base + q2 + i), pack_pair(R2));
 _mm_storeu_si128((__m128i*)(base + q3 + i), pack_pair(R3));
}

// j=0 special: v=0, u_lo=0、 u_hi=β_1 のみ非零。 1 PCLMUL のみ。
VPCLMUL inline void btf_fft_radix4_pair_j0(u64* base, int q, int q2, int q3, int i, __m256i u_hi_vec) {
 __m256i Q0= expand_pair(base + i);
 __m256i Q1= expand_pair(base + q + i);
 __m256i Q2= expand_pair(base + q2 + i);
 __m256i Q3= expand_pair(base + q3 + i);
 // Step 1 (v=0): A_0 = Q_0, A_1 = Q_1; A_2 = Q_0 ^ Q_2, A_3 = Q_1 ^ Q_3
 __m256i A2= _mm256_xor_si256(Q0, Q2);
 __m256i A3= _mm256_xor_si256(Q1, Q3);
 // Step 2 lower (u_lo=0): R_C00 = Q_0, R_C01 = Q_0 ^ Q_1
 __m256i R1= _mm256_xor_si256(Q0, Q1);
 // Step 2 upper: R_C10 = A_2 ^ u_hi A_3; R_C11 = R_C10 ^ A_3
 __m256i R2= _mm256_xor_si256(A2, clmul_reduce_pair(A3, u_hi_vec));
 __m256i R3= _mm256_xor_si256(R2, A3);
 _mm_storeu_si128((__m128i*)(base + i), pack_pair(Q0));
 _mm_storeu_si128((__m128i*)(base + q + i), pack_pair(R1));
 _mm_storeu_si128((__m128i*)(base + q2 + i), pack_pair(R2));
 _mm_storeu_si128((__m128i*)(base + q3 + i), pack_pair(R3));
}

// scalar fallback (q=1、 v=0 / u_lo=0 で skip 効果使うので残す)
inline void btf_fft_radix4_scalar(u64* base, int q, int q2, int q3, u64 v, u64 u_lo, u64 u_hi) {
 u64 Q0= base[0], Q1= base[q], Q2= base[q2], Q3= base[q3];
 u64 A0, A1, A2, A3;
 if(v == 0) {
  A0= Q0;
  A1= Q1;
 } else {
  A0= Q0 ^ mul(Q2, v);
  A1= Q1 ^ mul(Q3, v);
 }
 A2= A0 ^ Q2;
 A3= A1 ^ Q3;
 u64 R0, R1;
 if(u_lo == 0) R0= A0;
 else R0= A0 ^ mul(A1, u_lo);
 R1= R0 ^ A1;
 u64 R2= A2 ^ mul(A3, u_hi);
 u64 R3= R2 ^ A3;
 base[0]= R0;
 base[q]= R1;
 base[q2]= R2;
 base[q3]= R3;
}

// q=1 の SIMD 化 (v ≠ 0 かつ u_lo ≠ 0 を仮定; j ≥ 1 で常に成立)。
//   Step 1: pack (Q2, Q3) + broadcast v → 1 VPCLMULQDQ で (Q2·v, Q3·v)
//   Step 2: pack (A1, A3) + pack (u_lo, u_hi) → 1 VPCLMULQDQ で (A1·u_lo, A3·u_hi)
// q=1 専用なので q, q2, q3 引数は省略 (= 1, 2, 3)。
VPCLMUL inline void btf_fft_radix4_q1_simd(u64* base, __m256i v_vec, __m256i u_pair_vec) {
 u64 Q0= base[0], Q1= base[1], Q2= base[2], Q3= base[3];
 // Step 1: (Q2, Q3) × v、 Q2/Q3 はメモリ隣接なので expand_pair で 1 load + 1 permute
 __m256i Q23_vec= expand_pair(base + 2);
 __m256i p1= clmul_reduce_pair(Q23_vec, v_vec);
 u64 p_Q2v= (u64)_mm256_extract_epi64(p1, 0);
 u64 p_Q3v= (u64)_mm256_extract_epi64(p1, 2);
 u64 A0= Q0 ^ p_Q2v;
 u64 A1= Q1 ^ p_Q3v;
 u64 A2= A0 ^ Q2;
 u64 A3= A1 ^ Q3;
 // Step 2: (A1·u_lo, A3·u_hi)
 __m256i A13_vec= _mm256_set_epi64x(0, A3, 0, A1);
 __m256i p2= clmul_reduce_pair(A13_vec, u_pair_vec);
 u64 p_A1ulo= (u64)_mm256_extract_epi64(p2, 0);
 u64 p_A3uhi= (u64)_mm256_extract_epi64(p2, 2);
 u64 R0= A0 ^ p_A1ulo;
 u64 R1= R0 ^ A1;
 u64 R2= A2 ^ p_A3uhi;
 u64 R3= R2 ^ A3;
 base[0]= R0;
 base[1]= R1;
 base[2]= R2;
 base[3]= R3;
}

VPCLMUL inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 bc_to_lch(f.data(), n);
 const u64* M= nim_data.master.data();
 int k= d;
 while(k >= 2) {
  int unit= 1 << k;
  int q= unit >> 2;
  int q2= q << 1;
  int q3= q + q2;
  int num_blocks= n / unit;
  // j=0: v=0, u_lo=0、 u_hi=M[1]=β_1
  {
   u64* base= f.data();
   u64 u_hi= M[1];
   if(q == 1) {
    btf_fft_radix4_scalar(base, q, q2, q3, 0, 0, u_hi);
   } else {
    __m256i u_hi_vec= broadcast_ska(u_hi);
    int i= 0;
    for(; i + 1 < q; i+= 2) btf_fft_radix4_pair_j0(base, q, q2, q3, i, u_hi_vec);
    for(; i < q; ++i) btf_fft_radix4_scalar(base + i, q, q2, q3, 0, 0, u_hi);
   }
  }
  if(q == 1) {
   // q=1 専用: SIMD scalar (j > 0 で v ≠ 0, u_lo ≠ 0 を仮定)
   for(int j= 1; j < num_blocks; ++j) {
    u64* base= f.data() + (size_t)j * unit;
    __m256i v_vec= broadcast_ska(M[j]);
    __m256i u_pair_vec= pack_u_pair(M[2 * j], M[2 * j + 1]);
    btf_fft_radix4_q1_simd(base, v_vec, u_pair_vec);
   }
  } else {
   for(int j= 1; j < num_blocks; ++j) {
    u64* base= f.data() + (size_t)j * unit;
    u64 v= M[j];
    u64 u_lo= M[2 * j];
    u64 u_hi= M[2 * j + 1];
    __m256i v_vec= broadcast_ska(v);
    __m256i u_lo_vec= broadcast_ska(u_lo);
    __m256i u_hi_vec= broadcast_ska(u_hi);
    int i= 0;
    for(; i + 1 < q; i+= 2) btf_fft_radix4_pair(base, q, q2, q3, i, v_vec, u_lo_vec, u_hi_vec);
    for(; i < q; ++i) btf_fft_radix4_scalar(base + i, q, q2, q3, v, u_lo, u_hi);
   }
  }
  k-= 2;
 }
 // d 奇数なら最後に level 1 radix-2
 if(k == 1) {
  // level 1 (half=1) は SIMD 立たないので scalar
  int num_blocks= n >> 1;
  // j=0: ska=0
  {
   u64* base= f.data();
   base[1]^= base[0];
  }
  for(int j= 1; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * 2;
   u64 ska= M[j];
   u64 r= mul(base[1], ska);
   base[0]^= r;
   base[1]^= base[0];
  }
 }
}

// inverse butterfly group (general case)
VPCLMUL inline void btf_ifft_radix4_pair(u64* base, int q, int q2, int q3, int i, __m256i v_vec, __m256i u_lo_vec, __m256i u_hi_vec) {
 __m256i R0= expand_pair(base + i);
 __m256i R1= expand_pair(base + q + i);
 __m256i R2= expand_pair(base + q2 + i);
 __m256i R3= expand_pair(base + q3 + i);
 // Step 2 inverse upper: A_3 = R2 ^ R3, A_2 = R2 ^ u_hi A_3
 __m256i A3= _mm256_xor_si256(R2, R3);
 __m256i A2= _mm256_xor_si256(R2, clmul_reduce_pair(A3, u_hi_vec));
 // Step 2 inverse lower: A_1 = R0 ^ R1, A_0 = R0 ^ u_lo A_1
 __m256i A1= _mm256_xor_si256(R0, R1);
 __m256i A0= _mm256_xor_si256(R0, clmul_reduce_pair(A1, u_lo_vec));
 // Step 1 inverse: Q_2 = A_0 ^ A_2, Q_0 = A_0 ^ v Q_2; Q_3 = A_1 ^ A_3, Q_1 = A_1 ^ v Q_3
 __m256i Q2= _mm256_xor_si256(A0, A2);
 __m256i Q3= _mm256_xor_si256(A1, A3);
 __m256i Q0= _mm256_xor_si256(A0, clmul_reduce_pair(Q2, v_vec));
 __m256i Q1= _mm256_xor_si256(A1, clmul_reduce_pair(Q3, v_vec));
 _mm_storeu_si128((__m128i*)(base + i), pack_pair(Q0));
 _mm_storeu_si128((__m128i*)(base + q + i), pack_pair(Q1));
 _mm_storeu_si128((__m128i*)(base + q2 + i), pack_pair(Q2));
 _mm_storeu_si128((__m128i*)(base + q3 + i), pack_pair(Q3));
}

VPCLMUL inline void btf_ifft_radix4_pair_j0(u64* base, int q, int q2, int q3, int i, __m256i u_hi_vec) {
 __m256i R0= expand_pair(base + i);
 __m256i R1= expand_pair(base + q + i);
 __m256i R2= expand_pair(base + q2 + i);
 __m256i R3= expand_pair(base + q3 + i);
 // Step 2 inverse upper: A_3 = R2 ^ R3, A_2 = R2 ^ u_hi A_3
 __m256i A3= _mm256_xor_si256(R2, R3);
 __m256i A2= _mm256_xor_si256(R2, clmul_reduce_pair(A3, u_hi_vec));
 // Step 2 inverse lower (u_lo=0): A_1 = R0 ^ R1, A_0 = R0
 __m256i A1= _mm256_xor_si256(R0, R1);
 // Step 1 inverse (v=0): Q_2 = A_0 ^ A_2 = R_0 ^ A_2, Q_0 = A_0 = R_0; Q_3 = A_1 ^ A_3, Q_1 = A_1
 __m256i Q2= _mm256_xor_si256(R0, A2);
 __m256i Q3= _mm256_xor_si256(A1, A3);
 _mm_storeu_si128((__m128i*)(base + i), pack_pair(R0));
 _mm_storeu_si128((__m128i*)(base + q + i), pack_pair(A1));
 _mm_storeu_si128((__m128i*)(base + q2 + i), pack_pair(Q2));
 _mm_storeu_si128((__m128i*)(base + q3 + i), pack_pair(Q3));
}

inline void btf_ifft_radix4_scalar(u64* base, int q, int q2, int q3, u64 v, u64 u_lo, u64 u_hi) {
 u64 R0= base[0], R1= base[q], R2= base[q2], R3= base[q3];
 u64 A3= R2 ^ R3;
 u64 A2= R2 ^ mul(A3, u_hi);
 u64 A1= R0 ^ R1;
 u64 A0;
 if(u_lo == 0) A0= R0;
 else A0= R0 ^ mul(A1, u_lo);
 u64 Q2= A0 ^ A2;
 u64 Q3= A1 ^ A3;
 u64 Q0, Q1;
 if(v == 0) {
  Q0= A0;
  Q1= A1;
 } else {
  Q0= A0 ^ mul(Q2, v);
  Q1= A1 ^ mul(Q3, v);
 }
 base[0]= Q0;
 base[q]= Q1;
 base[q2]= Q2;
 base[q3]= Q3;
}

// q=1 SIMD inverse (j ≥ 1 で v ≠ 0, u_lo ≠ 0 を仮定)
//   Step 2 inverse: A_3 = R2^R3, A_2 = R2 ^ u_hi·A_3; A_1 = R0^R1, A_0 = R0 ^ u_lo·A_1
//     pack (A_3, A_1) + pack (u_hi, u_lo) → 1 VPCLMULQDQ
//   Step 1 inverse: Q_2 = A_0^A_2, Q_3 = A_1^A_3, Q_0 = A_0 ^ v·Q_2, Q_1 = A_1 ^ v·Q_3
//     pack (Q_2, Q_3) + broadcast v → 1 VPCLMULQDQ
VPCLMUL inline void btf_ifft_radix4_q1_simd(u64* base, __m256i v_vec, __m256i u_pair_inv_vec) {
 u64 R0= base[0], R1= base[1], R2= base[2], R3= base[3];
 // (A_3 = R2^R3 と A_1 = R0^R1 は scalar、 Step 2 mul は scalar pack 必要)
 u64 A3= R2 ^ R3;
 u64 A1= R0 ^ R1;
 // Step 2: (A3·u_hi, A1·u_lo)、 lane 0 = A3, lane 2 = A1。 u_pair_inv_vec = (u_hi, u_lo) 形式
 __m256i A31_vec= _mm256_set_epi64x(0, A1, 0, A3);
 __m256i p2= clmul_reduce_pair(A31_vec, u_pair_inv_vec);
 u64 p_A3uhi= (u64)_mm256_extract_epi64(p2, 0);
 u64 p_A1ulo= (u64)_mm256_extract_epi64(p2, 2);
 u64 A2= R2 ^ p_A3uhi;
 u64 A0= R0 ^ p_A1ulo;
 u64 Q2= A0 ^ A2;
 u64 Q3= A1 ^ A3;
 // Step 1: (Q2·v, Q3·v) with broadcast v
 __m256i Q23_vec= _mm256_set_epi64x(0, Q3, 0, Q2);
 __m256i p1= clmul_reduce_pair(Q23_vec, v_vec);
 u64 p_Q2v= (u64)_mm256_extract_epi64(p1, 0);
 u64 p_Q3v= (u64)_mm256_extract_epi64(p1, 2);
 u64 Q0= A0 ^ p_Q2v;
 u64 Q1= A1 ^ p_Q3v;
 base[0]= Q0;
 base[1]= Q1;
 base[2]= Q2;
 base[3]= Q3;
}

VPCLMUL inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 const u64* M= nim_data.master.data();
 int k_start;
 if(d & 1) {
  // level 1 radix-2 inverse 先行
  int num_blocks= n >> 1;
  {
   u64* base= f.data();
   base[1]^= base[0];
  }
  for(int j= 1; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * 2;
   u64 ska= M[j];
   base[1]^= base[0];
   base[0]^= mul(base[1], ska);
  }
  k_start= 3;
 } else {
  k_start= 2;
 }
 for(int level_hi= k_start; level_hi <= d; level_hi+= 2) {
  int unit= 1 << level_hi;
  int q= unit >> 2;
  int q2= q << 1;
  int q3= q + q2;
  int num_blocks= n / unit;
  {
   u64* base= f.data();
   u64 u_hi= M[1];
   if(q == 1) {
    btf_ifft_radix4_scalar(base, q, q2, q3, 0, 0, u_hi);
   } else {
    __m256i u_hi_vec= broadcast_ska(u_hi);
    int i= 0;
    for(; i + 1 < q; i+= 2) btf_ifft_radix4_pair_j0(base, q, q2, q3, i, u_hi_vec);
    for(; i < q; ++i) btf_ifft_radix4_scalar(base + i, q, q2, q3, 0, 0, u_hi);
   }
  }
  if(q == 1) {
   for(int j= 1; j < num_blocks; ++j) {
    u64* base= f.data() + (size_t)j * unit;
    __m256i v_vec= broadcast_ska(M[j]);
    // inverse: lane 0 = u_hi, lane 2 = u_lo (forward と入れ替わり)
    __m256i u_pair_inv_vec= pack_u_pair(M[2 * j + 1], M[2 * j]);
    btf_ifft_radix4_q1_simd(base, v_vec, u_pair_inv_vec);
   }
  } else {
   for(int j= 1; j < num_blocks; ++j) {
    u64* base= f.data() + (size_t)j * unit;
    u64 v= M[j];
    u64 u_lo= M[2 * j];
    u64 u_hi= M[2 * j + 1];
    __m256i v_vec= broadcast_ska(v);
    __m256i u_lo_vec= broadcast_ska(u_lo);
    __m256i u_hi_vec= broadcast_ska(u_hi);
    int i= 0;
    for(; i + 1 < q; i+= 2) btf_ifft_radix4_pair(base, q, q2, q3, i, v_vec, u_lo_vec, u_hi_vec);
    for(; i < q; ++i) btf_ifft_radix4_scalar(base + i, q, q2, q3, v, u_lo, u_hi);
   }
  }
 }
 bc_to_mono(f.data(), n);
}

VPCLMUL inline std::vector<u64> nim_convolution(std::vector<u64> f, std::vector<u64> g) {
 int n= (int)f.size(), m= (int)g.size();
 int s= (int)ceil_pow2(u32(n + m - 1));
 f.resize(s);
 g.resize(s);
 nim_data.init(msb(s));
 nim_fft(f);
 nim_fft(g);
 for(int i= 0; i < s; ++i) f[i]= mul(f[i], g[i]);
 nim_ifft(f);
 f.resize(n + m - 1);
 return f;
}

}  // namespace conv_f2_64_full_dif_radix4_v2

struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_full_dif_radix4_v2;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
