#pragma once
// GFNI (`_mm256_gf2p8affine_epi64_epi8`) ベースの squaring:
//   GFNI は per-byte の 8x8 F_2 affine 変換を 1 命令で実行可能。
//   sq の spread (bit i → bit 2i) は byte 局所な変換なので GFNI 1 発で full spread が出る。
//
// 構造:
//   1. v = __m256i に [a:a] を複製 (lane 0/1 とも入力 a)
//   2. M = [A_lo:A_hi] を行列引数に
//   3. _mm256_gf2p8affine_epi64_epi8(v, M, 0)
//      → lane 0: A_lo を適用 (各 byte の low nibble を 0,2,4,6 bit に spread)
//      → lane 1: A_hi を適用 (各 byte の high nibble を 0,2,4,6 bit に spread)
//   4. lane 0/1 を unpacklo で interleave → 128-bit 平方
//   5. reduction は scalar shift+XOR
//
// 必要な拡張: GFNI + AVX2 (Intel Ice Lake / AMD Zen4 以降)。
// AMD Zen3 までは GFNI 非対応 → 実機で動かない。SIMDe が ARM 等でエミュレーション可。
//
// PSHUFB 版との比較:
//   pshufb: 2 PSHUFB + 2 AND + 1 SRLI ≈ 3 cycle TPT
//   gfni:   1 gf2p8affine ≈ 1 cycle TPT, ~5 cycle latency
//   spread phase で ~3x throughput、reduction 込みで ~25% 改善見込み。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("gfni,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#ifdef USE_SIMDE
#include <simde/x86/gfni.h>
#endif
namespace gf2_64_sq_gfni {
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// 行列エンコード (Intel GFNI 仕様: A.byte[7-k] が出力 bit k に寄与する入力 bit のマスク):
//   A_LO: 入力 bit j (j∈[0,4)) → 出力 bit 2j (奇数 bit はゼロ)
//   A_HI: 入力 bit j (j∈[4,8)) → 出力 bit 2(j-4)
constexpr u64 A_LO= 0x0100020004000800ull;
constexpr u64 A_HI= 0x1000200040008000ull;
inline u64 sq(u64 a) {
 __m256i v= _mm256_set_epi64x(0, 0, (i64)a, (i64)a);
 __m256i M= _mm256_set_epi64x(0, 0, (i64)A_HI, (i64)A_LO);
 __m256i res= _mm256_gf2p8affine_epi64_epi8(v, M, 0);
 __m128i first128= _mm256_castsi256_si128(res);
 __m128i lo_s= first128;                       // bytes 0..7 (low_spread of each input byte)
 __m128i hi_s= _mm_srli_si128(first128, 8);    // bytes 8..15 (high_spread)
 __m128i squared= _mm_unpacklo_epi8(lo_s, hi_s);
 u64 lo= (u64)_mm_cvtsi128_si64(squared);
 u64 h= (u64)_mm_cvtsi128_si64(_mm_srli_si128(squared, 8));
 u64 d= h ^ (h << 1);
 return lo ^ RED[h >> 60] ^ d ^ (d << 3);
}
}  // namespace gf2_64_sq_gfni
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_gfni::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
