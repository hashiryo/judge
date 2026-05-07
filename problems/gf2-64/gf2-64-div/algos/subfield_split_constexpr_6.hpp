// subfield_split_constexpr_5 + reduction も __m256i で並列化:
//
// 5 では VPCLMULQDQ で 2 clmul を並列にしたが reduction はスカラー 2 回。
// 本版は reduction の主要部 (d = h^(h<<1), red1 = d^(d<<3), partial = lo^red1) を
// __m256i 演算で 2 lane 同時に実行。 RED[h>>60] の小さな byte 補正だけスカラーに残す。
//
// prod 配置 (256-bit):
//   要素 0 (bits 0..63)     = clmul(a32, N16) の low  → lo0
//   要素 1 (bits 64..127)   = clmul(a32, N16) の high → h0
//   要素 2 (bits 128..191)  = clmul(N,   N16) の low  → lo1
//   要素 3 (bits 192..255)  = clmul(N,   N16) の high → h1
//
// vectorize 戦略:
//   d_full   = prod ^ (prod << 1 epi64)         // 各 64-bit 要素ごとに左 1 シフト + XOR
//   red1     = d_full ^ (d_full << 3 epi64)
//   要素 1 = d0^(d0<<3), 要素 3 = d1^(d1<<3) (要素 0, 2 はゴミだが下流で使わない)
//   red1 を 128-bit lane 内で 8 byte 右シフトして要素 1,3 を low に持って来る
//   prod XOR red1_shifted で各 lane low に lo_k^red1_k が入る
//   最後に RED[h_k>>60] を scalar 加算
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
PCLMUL inline u64 inv(u64 a) {
 assert(a != 0);
 u64 a32= frob32(a);
 u64 N= mul(a, a32);
 u64 N16= frob16(N);
 // VPCLMULQDQ で 2 並列 clmul
 __m256i a_vec= _mm256_set_epi64x(0, N, 0, a32);
 __m256i b_vec= _mm256_set1_epi64x(N16);
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 // reduction の主要部を __m256i で並列実行
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 // 要素 1, 3 (= 各 lane の high 64 = d_k^(d_k<<3)) を各 lane の low 64 にシフト
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 // prod (各 lane low 64 = lo_k) と XOR → 各 lane low 64 = lo_k ^ red1_k
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 __m128i p0= _mm256_castsi256_si128(partial);
 __m128i p1= _mm256_extracti128_si256(partial, 1);
 // RED[h_k >> 60] を scalar で加算
 u64 g= (u64)p0[0] ^ RED[(u64)p0[1] >> 60];
 u64 NN16= (u64)p1[0] ^ RED[(u64)p1[1] >> 60];
 u64 b= embed_idx(INV_LOW[u16(NN16)]);
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
