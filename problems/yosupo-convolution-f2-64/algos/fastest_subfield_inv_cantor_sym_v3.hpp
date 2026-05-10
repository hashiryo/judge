#pragma once
#include "_common.hpp"
// cantor_sym_v2 派生: 余り 1 個を tail (i = half-1) ではなく **先頭 i = 1** に置く。
//
//   half=1 (len=2):                 i=0 のみ (PCLMUL なし)
//   half=2 (len=4):                 i=0 + i=1 scalar
//   half=N≥4:                       i=0 + i=1 scalar + pair(2,3)..pair(N-2,N-1)
//
// 利点として期待:
//   - i=0 (XOR-only) と i=1 (scalar PCLMUL, lat ~5) が dispatch 直後に発行され、
//     後続 pair ループとは独立に走るので OoO で latency を埋めやすい
//   - ブロック末尾で VPCLMULQDQ (mul2) で締めることで、 次ブロックの prefetch /
//     L1 fetch とパイプが整列しやすい (tail 版だと末尾 scalar PCLMUL のあと次の
//     ブロックの mul2 へ port 5 を切り替える必要があった)
//
// 命令数は v2 (tail 版) と同じ (port 5 トータル N/2)。 dispatch アラインメントだけ
// の比較で、 効果は微小と予想。
//
// Cantor 基底の正規化済み twiddle a[i] には恒等性
//     a[i][j + 2^(i-1)] = a[i][j] ^ 1
// が成り立つ (最終 basis 要素が 1 = s_W(β) = 1 の魔法)。 これを g0 = a[i] の下半分,
// g1 = a[i] の上半分 と置くと g1[k] = g0[k] ^ 1。
// よって fft の butterfly で必要な 2 つの mul は次のように共有可能:
//     hi * g0[k]              (= p)
//     hi * g1[k] = hi * g0[k] ^ hi (∵ GF(2) 上の分配律)
// すなわち 1 PCLMUL + 1 XOR で 2 butterfly 出力が両方計算できる。
// (元: 1 butterfly あたり 2 PCLMUL → 新: 1 butterfly あたり 1 PCLMUL)
//
// さらにこの「1 PCLMUL / butterfly」を 2 個まとめて VPCLMULQDQ (mul2) 1 命令で発行
// すれば、 **ピーク 2 butterfly / cycle** (元の 4 倍スループット) が狙える。
//
// IFFT 側の butterfly:
//     f0_new = lo*g1 + hi*g0 = lo*(g0+1) + hi*g0 = (lo+hi)*g0 + lo
//     f1_new = lo + hi
// → s = lo+hi をまず計算、 1 PCLMUL (s*g0) で済む。 同様に mul2 で 2 並列化。
//
// 各ブロック先頭の i=0 では g0[0]=0 (= 正規化済み twiddle 表の先頭は常に 0) なので
// PCLMUL も load も不要、 単に f0[0]=lo, f1[0]=lo^hi だけで済む。 これを
// cantor_sym のメインループから切り出し、 i=0 は単独処理 → 残り (i=1..half-1) を
// pair で潰す構造に変更。
//
//   half=1 (len=2):                 i=0 単独 (PCLMUL なし)
//   half=2 (len=4):                 i=0 単独 + i=1 単独 (1 scalar mul)
//   half=N≥4:                       i=0 単独 + pair(1,2)..pair(N-3,N-2) + i=N-1 単独
//
// ペアの mul2 は 2 lane 両方が「有効な PCLMUL」 になり (cantor_sym v1 では i=0 ペア
// で lane 0 が hi*0 = 0 の無駄演算だった) 命令効率が上がる。 ただし Ice Lake では
// VPCLMULQDQ も PCLMUL もスループット 1 / cycle (port 5 占有) なので命令数の差は
// 出ないが、 g0[0] の load 削減と分岐レス化でわずかに前段 IPC が改善する見込み。
//
// その他は raw_u64 と同一: subfield-split inv, constexpr CHAIN, Gray code init。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"
#include "../../gf2-64/_shared/frob.hpp"

namespace conv_f2_64_cantor_sym_v3 {

using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// constexpr GF(2^64) mul / sq (chain を .rodata に焼くため)
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
// Artin-Schreier 連鎖 c[k] = P^k(2), P(x) = x²+x。
// 実測で c[62] = 1, c[63] = 0 なので非零 basis は c[0..62] の 63 個。
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();

// =============================================================================
// subfield-split inv (constchain 版と同一)
// =============================================================================
constexpr auto INV_LOW= []() {
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {}, T_hi[256]= {};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j)
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
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
  for(int half= 0; half < 2; ++half)
   for(int i= 0; i < 256; ++i) {
    u64 v= 0;
    for(int b= 0; b < 8; ++b)
     if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
    t[half][i]= v;
   }
  return t;
 }();
 return EMBED[0][u8(idx)] ^ EMBED[1][idx >> 8];
}

inline const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);

VPCLMUL inline void mul2(const __m256i& a_vec, const __m256i& b_vec, u64& r0, u64& r1) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 __m256i result= _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
 r0= _mm256_extract_epi64(result, 0);
 r1= _mm256_extract_epi64(result, 2);
}

VPCLMUL inline u64 inv(u64 a) {
 assert(a != 0);
 u64 a32= frob32(a);
 u64 N= mul(a, a32);
 u64 N16= frob16(N);
 mul2(_mm256_set_epi64x(0, N, 0, a32), _mm256_set1_epi64x(N16), a, N16);
 return mul(embed_idx(INV_LOW[u16(N16)]), a);
}

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// Twiddle: u64 配列で直接持つ (gf2 ラッパなし)
// =============================================================================
struct nim_fft_data {
 std::vector<std::vector<u64>> a;
 VPCLMUL void init(int n) {
  if((int)a.size() > n) return;
  a.resize(n + 1);
  std::vector<u64> b(n);
  for(int k= 0; k < n; ++k) b[k]= CHAIN[CHAIN_LEN - n + k];
  for(int i= n;; --i) {
   if(i > 0) {
    u64 inv_b= inv(b.back());
    for(u64& x : b) x= mul(x, inv_b);
   }
   auto& na= a[i];
   na.resize(1 << i);
   na[0]= 0;
   for(int j= 1; j < (1 << i); ++j) {
    int k= __builtin_ctz(j);
    na[j]= na[j ^ (1 << k)] ^ b[k];
   }
   if(i == 0) break;
   b.pop_back();
   for(u64& x : b) x= sq(x) ^ x;
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// nim FFT (u64 そのまま、 + は ^、 * は mul / mul2)
// =============================================================================
VPCLMUL inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // Phase A: XOR-only butterfly + deinterleave
 while(len > 1) {
  for(int l= 0; l < n; l+= len) {
   for(int m= len / 4; m >= 1; m/= 2) {
    for(int t= 0; t < len; t+= m * 4) {
     for(int i= 0; i < m; ++i) {
      u64 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b ^ c ^ d;
      f[l + t + m * 2 + i]= c ^ d;
     }
    }
   }
   for(int i= 0; i < len / 2; ++i) {
    f2[l + i]= f[l + i * 2];
    f2[l + i + len / 2]= f[l + i * 2 + 1];
   }
  }
  std::swap(f, f2);
  len/= 2;
 }
 // Phase B: Cantor 対称性 + g0[0]=0 スキップ + i=1 を scalar、 残り (i=2..) を pair
 while(len < n) {
  len*= 2;
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g0= g.data();
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   // i = 0: g0[0] = 0 → PCLMUL なし
   {
    u64 lo= f0[0], hi= f1[0];
    f0[0]= lo;
    f1[0]= lo ^ hi;
   }
   // i = 1: scalar mul (パイプ整列のため pair より前に置く)
   if(half >= 2) {
    u64 lo= f0[1], hi= f1[1];
    u64 p= mul(hi, g0[1]);
    u64 t= lo ^ p;
    f0[1]= t;
    f1[1]= t ^ hi;
   }
   // i = 2..half-1 を pair で mul2
   for(int i= 2; i + 1 < half; i+= 2) {
    u64 lo0= f0[i], hi0= f1[i];
    u64 lo1= f0[i + 1], hi1= f1[i + 1];
    u64 p0, p1;
    mul2(_mm256_set_epi64x(0, hi1, 0, hi0), _mm256_set_epi64x(0, g0[i + 1], 0, g0[i]), p0, p1);
    u64 t0= lo0 ^ p0;
    u64 t1= lo1 ^ p1;
    f0[i]= t0;
    f1[i]= t0 ^ hi0;
    f0[i + 1]= t1;
    f1[i + 1]= t1 ^ hi1;
   }
  }
 }
}

VPCLMUL inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // IFFT: Cantor 対称性 + g0[0]=0 スキップ + i=1 を scalar、 残り (i=2..) を pair
 while(len > 1) {
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g0= g.data();
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   // i = 0: s*g0 = 0
   {
    u64 lo= f0[0], hi= f1[0];
    f0[0]= lo;
    f1[0]= lo ^ hi;
   }
   // i = 1: scalar mul
   if(half >= 2) {
    u64 lo= f0[1], hi= f1[1];
    u64 s= lo ^ hi;
    u64 p= mul(s, g0[1]);
    f0[1]= lo ^ p;
    f1[1]= s;
   }
   for(int i= 2; i + 1 < half; i+= 2) {
    u64 lo0= f0[i], hi0= f1[i];
    u64 lo1= f0[i + 1], hi1= f1[i + 1];
    u64 s0= lo0 ^ hi0;
    u64 s1= lo1 ^ hi1;
    u64 p0, p1;
    mul2(_mm256_set_epi64x(0, s1, 0, s0), _mm256_set_epi64x(0, g0[i + 1], 0, g0[i]), p0, p1);
    f0[i]= lo0 ^ p0;
    f0[i + 1]= lo1 ^ p1;
    f1[i]= s0;
    f1[i + 1]= s1;
   }
  }
  len/= 2;
 }
 while(len < n) {
  len*= 2;
  for(int l= 0; l < n; l+= len) {
   for(int i= 0; i < len / 2; ++i) {
    f2[l + i * 2]= f[l + i];
    f2[l + i * 2 + 1]= f[l + i + len / 2];
   }
  }
  std::swap(f, f2);
  for(int l= 0; l < n; l+= len) {
   for(int m= 1; m <= len / 4; m*= 2) {
    for(int t= 0; t < len; t+= m * 4) {
     for(int i= 0; i < m; ++i) {
      u64 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b ^ c;
      f[l + t + m * 2 + i]= c ^ d;
     }
    }
   }
  }
 }
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

}  // namespace conv_f2_64_cantor_sym_v3

struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_cantor_sym_v3;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
