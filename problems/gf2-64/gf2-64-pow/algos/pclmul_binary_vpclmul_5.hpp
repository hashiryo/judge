#pragma once
// pclmul_binary_vpclmul_4 + state を __m256i のまま iter 間で保持:
//
// vpclmul_4 では各 iter で res, a を scalar で持ち、 step 関数の中で
// __m256i に詰め直して、 終了時に scalar に戻していた。
// 本版は state __m256i を直接持ち回す:
//   state の lane 0 low 64 = res, lane 1 low 64 = a (high 64 はゴミ、 clmul imm=0 で無視)
//
// step iter (bit=1):
//   1. b_vec = permute4x64 で element 2 (= a) を全 lane に broadcast
//   2. clmul(state, b_vec, 0)
//   3. reduction (vpclmul_4 と同じ shift+xor + PSHUFB RED)
//   4. state = final_vec (next iter の a_vec として直接再利用可能なレイアウト)
//
// bit=0 iter:
//   a を一旦 scalar に抽出 → scalar sq → element 2 に insert で書き戻し
//   (vpclmul で sq だけ計算するのは無駄)
//
// 利点: bit=1 の step iter で a_vec/b_vec の scalar→vector 詰め直しコスト (各 ~5 cycle) を削減。
//        connecting iter 間で extract/repack が不要。
// 欠点: bit=0 で extract/insert オーバーヘッド (~5 cycle) が乗る。
//        random e (popcount ~32) では概ね break-even、 dense_e で有利、 sparse e で不利。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary_vpclmul_5 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
// state: lane 0 low 64 = res, lane 1 low 64 = a, high 64 はゴミ (clmul imm=0 で無視)
VPCLMUL inline __m256i step_vec(__m256i state) {
 // b_vec: element 2 (= a) を 全 lane に broadcast
 __m256i b_vec= _mm256_permute4x64_epi64(state, _MM_SHUFFLE(2, 2, 2, 2));
 __m256i prod= _mm256_clmulepi64_epi128(state, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 __m256i h_idx= _mm256_srli_epi64(partial, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 return _mm256_xor_si256(partial, red_vec);
}
inline u64 pow_binary(u64 a, u64 e) {
 if(!e) return 1;
 // 初期 state: res=1 (element 0), a (element 2)
 __m256i state= _mm256_set_epi64x(0, (i64)a, 0, 1);
 while(e) {
  if(e & 1) {
   if(e == 1) {
    // 最終 iter: res ← mul(res, a) のみ必要 (sq(a) は捨てる)。 step_vec で両方計算して res だけ取る。
    state= step_vec(state);
    break;
   }
   state= step_vec(state);
  } else {
   // a を scalar 取り出し → sq → element 2 に書き戻し
   u64 a_scalar= (u64)_mm256_extracti128_si256(state, 1)[0];
   a_scalar= sq(a_scalar);
   state= _mm256_insert_epi64(state, (i64)a_scalar, 2);
  }
  e>>= 1;
 }
 return (u64)_mm256_castsi256_si128(state)[0];
}
}  // namespace gf2_64_pow_pclmul_binary_vpclmul_5
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary_vpclmul_5::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
