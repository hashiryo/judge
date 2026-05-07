#pragma once
// gfni_linmap.hpp の改良版 (column-major + accumulator):
//
// 元 gfni_linmap.hpp は p (input byte) を外側にしたため、 lane 4 個から 1 byte ずつ
// 抽出して u64 を組む処理を p ごとに繰り返し → memory store + scalar XOR が支配的。
//
// 本版は J (input column block) を外側に、 4 row block を 1 __m256i にまとめてロード:
//   y_I = ⊕_J M_{I,J} · x_J   (output byte I, J=0..7 で集約)
// J ループ毎に
//   - x_J を YMM 全 byte に broadcast
//   - 4 row block (M_{0..3,J} or M_{4..7,J}) を 1 __m256i load
//   - gf2p8affine で 4 lane 同時に M_pj · x_J 計算
//   - YMM XOR で accumulator に蓄積
// 蓄積後、 各 lane は [output_byte_I × 8] (broadcast) になる → 1 shuffle で gather。
//
// per call で memory store + scalar 抽出が不要 → gfni_linmap.hpp より高速見込み。
//
// アイデア出典: 別 LLM の提案。 ただし行列 encoding が逆順だったので本実装で修正。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("gfni,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#ifdef USE_SIMDE
#include <simde/x86/gfni.h>
#endif
namespace gf2_64_sq_gfni_linmap_accum {
constexpr u64 sq_naive(u64 a) {
 u64 lo= 0, hi= 0;
 for(int i= 0; i < 64; ++i)
  if((a >> i) & 1) {
   int j= 2 * i;
   if(j < 64) lo^= u64(1) << j;
   else hi^= u64(1) << (j - 64);
  }
 for(int i= 63; i >= 0; --i)
  if((hi >> i) & 1) {
   hi^= u64(1) << i;
   lo^= IRRED_LOW << i;
   if(i > 0) hi^= IRRED_LOW >> (64 - i);
  }
 return lo;
}
// 8x8 サブ行列 M_{I,J} を canonical gf2p8affine encoding で u64 へ
//   byte[7-k] = row k of M_{I,J} (= bit q が M_{I,J}[k][q])
constexpr u64 encode_submat(int I, int J) {
 u64 enc= 0;
 for(int k= 0; k < 8; ++k) {
  u8 row_bits= 0;
  for(int q= 0; q < 8; ++q) {
   u64 sq_val= sq_naive(u64(1) << (8 * J + q));
   if((sq_val >> (8 * I + k)) & 1) row_bits|= u8(1) << q;
  }
  enc|= u64(row_bits) << (8 * (7 - k));
 }
 return enc;
}
// BLOCKS[J*8 + I] = サブ行列 (column-major: J 外、I 内)
struct alignas(64) BlockTable {
 u64 b[64];
};
constexpr BlockTable BLOCKS= []() {
 BlockTable t{};
 for(int I= 0; I < 8; ++I)
  for(int J= 0; J < 8; ++J) t.b[J * 8 + I]= encode_submat(I, J);
 return t;
}();
inline u64 sq(u64 a) {
 // x_bcast: 全 lane が a (= 各 lane の bytes 0..7 が a の byte 0..7)
 __m256i x_bcast= _mm256_set1_epi64x((i64)a);
 __m256i y_lo= _mm256_setzero_si256();  // I=0..3 の累積
 __m256i y_hi= _mm256_setzero_si256();  // I=4..7 の累積
 // PSHUFB で a の byte J を全 32 byte に broadcast するためのマスク
 // 各 J ∈ [0, 8) について「すべての byte が J」の 32 byte mask
 alignas(32) static const u8 BCAST_MASKS[8][32]= {
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4},
  {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5},
  {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
  {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
 };
#pragma GCC unroll 8
 for(int J= 0; J < 8; ++J) {
  __m256i x_J= _mm256_shuffle_epi8(x_bcast, _mm256_load_si256((const __m256i*)BCAST_MASKS[J]));
  __m256i M_lo= _mm256_load_si256((const __m256i*)&BLOCKS.b[J * 8]);     // M_{0..3, J}
  __m256i M_hi= _mm256_load_si256((const __m256i*)&BLOCKS.b[J * 8 + 4]); // M_{4..7, J}
  __m256i p_lo= _mm256_gf2p8affine_epi64_epi8(x_J, M_lo, 0);
  __m256i p_hi= _mm256_gf2p8affine_epi64_epi8(x_J, M_hi, 0);
  y_lo= _mm256_xor_si256(y_lo, p_lo);
  y_hi= _mm256_xor_si256(y_hi, p_hi);
 }
 // 各 lane は [output_byte_I × 8] になっている → byte 0 (と byte 8 = lane 1 byte 0) を gather
 alignas(32) static const u8 GATHER[32]= {
  0, 8, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
  0, 8, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
 };
 __m256i mask= _mm256_load_si256((const __m256i*)GATHER);
 __m256i y_lo_g= _mm256_shuffle_epi8(y_lo, mask);
 __m256i y_hi_g= _mm256_shuffle_epi8(y_hi, mask);
 // 下位 128 bit に lanes 0,1 / 上位 128 bit に lanes 2,3 の byte 0 が揃う
 u64 b01= (u64)_mm_cvtsi128_si64(_mm256_castsi256_si128(y_lo_g)) & 0xFFFF;
 u64 b23= (u64)_mm_cvtsi128_si64(_mm256_extracti128_si256(y_lo_g, 1)) & 0xFFFF;
 u64 b45= (u64)_mm_cvtsi128_si64(_mm256_castsi256_si128(y_hi_g)) & 0xFFFF;
 u64 b67= (u64)_mm_cvtsi128_si64(_mm256_extracti128_si256(y_hi_g, 1)) & 0xFFFF;
 return b01 | (b23 << 16) | (b45 << 32) | (b67 << 48);
}
}  // namespace gf2_64_sq_gfni_linmap_accum
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_gfni_linmap_accum::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
