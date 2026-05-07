#pragma once
// frobenius_byte_2 と同じ「線形写像をテーブル経由で評価」を GFNI で実装したもの。
//
// 設計:
//   sq の 64x64 線形写像 M を 8x8 サブ行列 M_pj (p=入力 byte 位置, j=出力 byte 位置) に分解。
//     output_byte_j = ⊕_{p=0..7} M_pj · b_p
//   M_pj は 8x8 F_2 行列 → gf2p8affine の行列引数 (u64) にエンコードして compile time に保存。
//
// per call の使い方 (AVX2 GFNI = 4 lanes):
//   ある (p, j) ペアに対して: 入力 byte b_p を全 lane に broadcast、 lane に M_pj を置いて
//   gf2p8affine 適用すると lane 全 8 byte が [M_pj · b_p] (= 同じ値の 8 個 broadcast)。
//   1 call で 4 lane 分 = 4 個の (p, j) 寄与を同時計算。
//
// 全体: 8 input byte × 8 output byte = 64 (p, j) ペア。 4 lane / call で 16 calls 必要。
//   各 call から 1 byte ずつ抽出して output_byte に XOR、 最終的に u64 result を組み立てる。
//
// 性能予想: 16 GFNI calls + extraction overhead ≈ 50-60 cycle。
//   frobenius_byte_2 (8 lookups + 7 XOR ≈ 10 cycle) より明確に遅い。
//   GFNI が wide な線形写像に対して構造的に不向きなことを示すサンプル実装。
//   逆に言うと「線形写像をテーブル不要で表現できる」点が GFNI 版の魅力。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("gfni,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#ifdef USE_SIMDE
#include <simde/x86/gfni.h>
#endif
namespace gf2_64_sq_gfni_linmap {
// 素朴 sq (前計算用のみ、 ランタイムでは使わない)
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
// (p, j) サブ行列を gf2p8affine 仕様の u64 にエンコード:
//   行列 M_pj[k][q] = "入力 bit (8p+q) が出力 bit (8j+k) に寄与するか"
//   gf2p8affine: A.byte[7-k] が出力 bit k に対応する入力マスク (bit q が M_pj[k][q])
constexpr u64 encode_submatrix(int p, int j) {
 u64 enc= 0;
 for(int k= 0; k < 8; ++k) {
  u8 row_bits= 0;
  for(int q= 0; q < 8; ++q) {
   u64 sq_val= sq_naive(u64(1) << (8 * p + q));
   if((sq_val >> (8 * j + k)) & 1) row_bits|= u8(1) << q;
  }
  enc|= u64(row_bits) << (8 * (7 - k));
 }
 return enc;
}
// 64 個の 8x8 サブ行列を compile time に展開
constexpr auto MATS= []() {
 array<array<u64, 8>, 8> m{};
 for(int p= 0; p < 8; ++p)
  for(int j= 0; j < 8; ++j) m[p][j]= encode_submatrix(p, j);
 return m;
}();
inline u64 sq(u64 a) {
 u64 result= 0;
 // 各入力 byte p に対し、 出力 byte j=0..7 への寄与を 2 GFNI calls で計算
 for(int p= 0; p < 8; ++p) {
  u8 b_p= u8(a >> (8 * p));
  __m256i v= _mm256_set1_epi64x((i64)(u64(b_p) * 0x0101010101010101ull));
  // 4 lane に M_p,0..M_p,3 を、 別 call で M_p,4..M_p,7 を置く
  __m256i M0= _mm256_set_epi64x((i64)MATS[p][3], (i64)MATS[p][2], (i64)MATS[p][1], (i64)MATS[p][0]);
  __m256i M1= _mm256_set_epi64x((i64)MATS[p][7], (i64)MATS[p][6], (i64)MATS[p][5], (i64)MATS[p][4]);
  __m256i r0= _mm256_gf2p8affine_epi64_epi8(v, M0, 0);
  __m256i r1= _mm256_gf2p8affine_epi64_epi8(v, M1, 0);
  // 各 lane の 8 byte は同じ値が broadcast されている → byte 0 だけ取れば十分。
  // メモリ経由で 4 lane × u64 を取り出し、 各 u64 の最下位 byte が contribution。
  alignas(32) u64 buf0[4], buf1[4];
  _mm256_store_si256((__m256i*)buf0, r0);
  _mm256_store_si256((__m256i*)buf1, r1);
  u64 contribution= (u64(u8(buf0[0])) << 0) | (u64(u8(buf0[1])) << 8) | (u64(u8(buf0[2])) << 16) | (u64(u8(buf0[3])) << 24) | (u64(u8(buf1[0])) << 32) | (u64(u8(buf1[1])) << 40) | (u64(u8(buf1[2])) << 48) | (u64(u8(buf1[3])) << 56);
  result^= contribution;
 }
 return result;
}
}  // namespace gf2_64_sq_gfni_linmap
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_gfni_linmap::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
