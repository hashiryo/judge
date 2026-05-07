// subfield_split_constexpr_5: 独立な 2 mul を VPCLMULQDQ で並列実行
//
// 元の inv() 内で:
//   g    = mul(a32, N16)        ← clmul(a32, N16) + reduction
//   tmp  = mul(N,   N16)        ← clmul(N,   N16) + reduction (u16 のみ使用)
// は a32, N, N16 が出揃った後、互いに独立。
// VPCLMULQDQ (AVX2) の `_mm256_clmulepi64_epi128` で 2 PCLMUL を 1 命令に。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降)。
// AMD EPYC 7763 = Zen3 → VPCLMULQDQ 対応 (GFNI と違って dashboard で動く期待)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"

using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// === メイン: low16 → low16 の inverse テーブル ===
constexpr auto INV_LOW= []() {
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {};
 u16 T_hi[256]= {};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j) {
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
   }
  }
  T_lo[v]= lo;
  T_hi[v]= hi;
 }
 u16 nat[65535];
 for(u16 i= 0, cur= 1; i < 65535; ++i) nat[i]= cur, cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 array<u16, 65536> t{};
 u16 lo1= T_lo[1] ^ T_hi[0];
 t[lo1]= lo1;
 for(uint32_t k= 1; k <= 32767; ++k) {
  u16 nk= nat[k];
  u16 nik= nat[65535 - k];
  u16 lo_k= T_lo[u8(nk)] ^ T_hi[nk >> 8];
  u16 lo_ik= T_lo[u8(nik)] ^ T_hi[nik >> 8];
  t[lo_k]= lo_ik;
  t[lo_ik]= lo_k;
 }
 return t;
}();
// === embed テーブル(low16 → 64bit、SUBFIELD_BASIS による)===
constexpr inline u64 embed_idx(u16 idx) {
 static constexpr auto EMBED= []() {
  u64 SUBFIELD_BASIS[]= {1ULL, 6899425322512154626ULL, 12712641506861907972ULL, 12687683756412895240ULL, 13108774640850436112ULL, 1196746230653255712ULL, 13779846473293824064ULL, 1136705091741089920ULL, 13132935623751303424ULL, 12256911237861802496ULL, 1968662052679910400ULL, 13476734309037115392ULL, 31478309824172032ULL, 5397840376063860736ULL, 18145356609018085376ULL, 2133828226494464000ULL};
  array<array<u64, 256>, 2> t{};
  for(int half= 0; half < 2; ++half) {
   for(int i= 0; i < 256; ++i) {
    u64 v= 0;
    for(int b= 0; b < 8; ++b)
     if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
    t[half][i]= v;
   }
  }
  return t;
 }();
 return EMBED[0][u8(idx)] ^ EMBED[1][idx >> 8];
}
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
VPCLMUL inline u64 inv(u64 a) {
 assert(a != 0);
 u64 a32= frob32(a);
 u64 N= mul(a, a32);
 u64 N16= frob16(N);
 // VPCLMULQDQ で 2 並列 clmul:
 //   lane 0: clmul(a32, N16) → 後で reduce → g
 //   lane 1: clmul(N,   N16) → 後で reduce → mul(N, N16)
 // _mm256_set_epi64x(e3, e2, e1, e0): 128-bit lane 0 = [e1:e0], lane 1 = [e3:e2]
 // PCLMUL imm8=0 は各 lane の low64 を mul する → e0, e2 だけ意味あり。
 __m256i a_vec= _mm256_set_epi64x(0, N, 0, a32);
 __m256i b_vec= _mm256_set1_epi64x(N16);
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m128i p0= _mm256_castsi256_si128(prod);       // clmul(a32, N16)
 __m128i p1= _mm256_extracti128_si256(prod, 1);  // clmul(N, N16)
 // reduce each (mul.hpp と同形)
 u64 g_h= (u64)p0[1];
 u64 d_g= g_h ^ (g_h << 1);
 u64 g= u64(p0[0]) ^ RED[g_h >> 60] ^ d_g ^ (d_g << 3);
 u16 t_h= p1[1], d_t= t_h ^ (t_h << 1);
 u16 NN16= u16(p1[0]) ^ RED[p1[1] >> 60] ^ d_t ^ (d_t << 3);
 // 残り
 u64 b= embed_idx(INV_LOW[NN16]);
 return mul(b, g);
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
