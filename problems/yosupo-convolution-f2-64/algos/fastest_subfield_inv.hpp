#pragma once
#include "_common.hpp"
// fastest_vpclmul ベース + 自作 inv (gf2-64-div/subfield_split_constexpr_7) 移植版。
//
// 改良点:
//   1. gf2::inv を pow(2^64-2) (~128 mul/sq 連鎖) から
//      F_{2^16} 部分体に落として lookup 1 回で終わらせる subfield_split に置換。
//      → nim_fft_data::init の inv が level 数だけ呼ばれるので precompute が短縮。
//   2. gf2::mul は gf2_64_pclmul::mul (PCLMUL + scalar reduction、 同一 idiom) を使用。
//   3. Phase B (fft) と twiddle 段 (ifft) の butterfly mul を VPCLMULQDQ 1 命令 +
//      __m256i 並列 reduction (subfield_split_constexpr_7 の mul2 と同じ idiom) に
//      置換。 fastest_vpclmul の mul_pair ではスカラー reduction だった部分も SIMD 化。
//
// アルゴリズム自体 (additive nim FFT) は fastest_vpclmul.hpp と同一。
//   基底列 b[i+1] = b[i]^2 + b[i] (Artin-Schreier) で additive subspace を生成し
//   そこを評価点に取る O((n+m) log(n+m)) FFT + reduction polynomial
//   x^64 + x^4 + x^3 + x + 1 の上での乗算。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"
#include "../../gf2-64/_shared/frob.hpp"

namespace conv_f2_64_subfield_inv {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// subfield_split_constexpr_7 から inv 関連定数 / ヘルパを移植
// =============================================================================
// F_{2^16}^* の log table から作る "inv low" map: F_{2^16} 内で x ↔ x^{-1} を直引き。
// 添字は norm N16 = N^{2^16} の F_{2^16} 表現の下位 16-bit。
constexpr auto INV_LOW= []() {
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {};
 u16 T_hi[256]= {};
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

// F_{2^16} の元 idx を F_{2^64} に embed (256-entry × 2 byte の split LUT)。
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

// VPCLMULQDQ で 2 mul + 並列 reduction を 1 関数で実行。
// a_vec / b_vec は lane 0 (要素 0) と lane 1 (要素 2) に被乗数を配置 (要素 1, 3 は 0)。
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

// 高速 inv: a^{-1} = N^{-1} * a^{2^32}*a (norm 経路) を subfield 内 lookup で。
//   N    = a * a^{2^32}                (∈ F_{2^32})
//   N16  = N^{2^16}                    (∈ F_{2^16})
//   N^{-1} は F_{2^16}^* で逆元 lookup → embed_idx で F_{2^64} に持ち上げ
//   a^{-1} = (a * a32) を mul2 でまとめて N16 に掛けつつ最終 mul で完成。
VPCLMUL inline u64 inv(u64 a) {
 assert(a != 0);
 u64 a32= frob32(a);
 u64 N= mul(a, a32);
 u64 N16= frob16(N);
 mul2(_mm256_set_epi64x(0, N, 0, a32), _mm256_set1_epi64x(N16), a, N16);
 return mul(embed_idx(INV_LOW[u16(N16)]), a);
}

// =============================================================================
// gf2 wrapper: + は XOR, * は PCLMUL ベース mul, inv は subfield_split。
// =============================================================================
struct gf2 {
 u64 v;
 gf2(): v(0) {}
 gf2(u64 x): v(x) {}
 gf2 operator+(const gf2& r) const { return gf2(v ^ r.v); }
 gf2 operator-(const gf2& r) const { return gf2(v ^ r.v); }
 gf2& operator+=(const gf2& r) { v^= r.v; return *this; }
 gf2& operator-=(const gf2& r) { v^= r.v; return *this; }
 PCLMUL gf2& operator*=(const gf2& r) {
  v= mul(v, r.v);
  return *this;
 }
 PCLMUL gf2 operator*(const gf2& r) const { return gf2(mul(v, r.v)); }
 VPCLMUL gf2 inverse() const { return gf2(inv(v)); }
 bool operator==(const gf2& r) const { return v == r.v; }
};

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// Twiddle 表 (additive FFT 用): a[i] は level i の 2^i 個の twiddle。
//   inv 計算は新しい subfield-split inv を利用。
// =============================================================================
struct nim_fft_data {
 std::vector<std::vector<gf2>> a;
 VPCLMUL void init(int n) {
  if((int)a.size() > n) return;
  a.resize(n + 1);
  std::vector<gf2> b;
  b.push_back(gf2(2));
  while(true) {
   gf2 x= b.back();
   if(x.v == 0) break;
   b.push_back(x * x + x);  // Artin-Schreier
  }
  b.pop_back();
  b.erase(b.begin(), b.end() - n);
  for(int i= n;; --i) {
   std::vector<gf2>& na= a[i];
   na.resize(1 << i);
   gf2 inv_b= b.back().inverse();
   for(gf2& x : b) x*= inv_b;
   for(int j= 0; j < (1 << i); ++j) {
    for(int k= 0; k < i; ++k)
     if(j >> k & 1) na[j]+= b[k];
   }
   if(i == 0) break;
   b.pop_back();
   for(gf2& x : b) x= x * x + x;
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// nim FFT / IFFT (fastest_vpclmul のループ構造を踏襲、mul_pair → mul2 に置換)
// =============================================================================
VPCLMUL inline void nim_fft(std::vector<gf2>& f) {
 int n= (int)f.size();
 std::vector<gf2> f2(n);
 int len= n;
 // Phase A: 乗算なし。 butterfly + deinterleave の前処理段。
 while(len > 1) {
  for(int l= 0; l < n; l+= len) {
   for(int m= len / 4; m >= 1; m/= 2) {
    for(int t= 0; t < len; t+= m * 4) {
     for(int i= 0; i < m; ++i) {
      gf2 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b + c + d;
      f[l + t + m * 2 + i]= c + d;
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
 // Phase B: in-place + VPCLMULQDQ で (hi*g0, hi*g1) を 1 命令並列 mul + 並列 reduction
 while(len < n) {
  len*= 2;
  const std::vector<gf2>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const gf2* g0= g.data();
   const gf2* g1= g.data() + half;
   gf2* f0= f.data() + l;
   gf2* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i].v, hi= f1[i].v;
    u64 p0, p1;
    mul2(_mm256_set_epi64x(0, hi, 0, hi), _mm256_set_epi64x(0, g1[i].v, 0, g0[i].v), p0, p1);
    f0[i]= gf2(lo ^ p0);
    f1[i]= gf2(lo ^ p1);
   }
  }
 }
}

VPCLMUL inline void nim_ifft(std::vector<gf2>& f) {
 int n= (int)f.size();
 std::vector<gf2> f2(n);
 int len= n;
 // 逆 Phase B: (lo*g1, hi*g0) を 1 VPCLMULQDQ で並列 mul
 while(len > 1) {
  const std::vector<gf2>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const gf2* g0= g.data();
   const gf2* g1= g.data() + half;
   gf2* f0= f.data() + l;
   gf2* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i].v, hi= f1[i].v;
    u64 p0, p1;
    mul2(_mm256_set_epi64x(0, hi, 0, lo), _mm256_set_epi64x(0, g0[i].v, 0, g1[i].v), p0, p1);
    f0[i]= gf2(p0 ^ p1);
    f1[i]= gf2(lo ^ hi);
   }
  }
  len/= 2;
 }
 // 逆 Phase A: deinterleave + butterfly (乗算なし)
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
      gf2 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b + c;
      f[l + t + m * 2 + i]= c + d;
     }
    }
   }
  }
 }
}

VPCLMUL inline std::vector<gf2> nim_convolution(std::vector<gf2> f, std::vector<gf2> g) {
 int n= (int)f.size(), m= (int)g.size();
 int s= (int)ceil_pow2(u32(n + m - 1));
 f.resize(s);
 g.resize(s);
 nim_data.init(msb(s));
 nim_fft(f);
 nim_fft(g);
 for(int i= 0; i < s; ++i) f[i]*= g[i];
 nim_ifft(f);
 f.resize(n + m - 1);
 return f;
}

}  // namespace conv_f2_64_subfield_inv

struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_subfield_inv;
  std::vector<gf2> a(n), b(m);
  for(int i= 0; i < n; ++i) a[i]= gf2(a_in[i]);
  for(int i= 0; i < m; ++i) b[i]= gf2(b_in[i]);
  auto c= nim_convolution(std::move(a), std::move(b));
  std::vector<u64> out(c.size());
  for(size_t k= 0; k < c.size(); ++k) out[k]= c[k].v;
  return out;
 }
};
