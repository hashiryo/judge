#pragma once
#include "_common.hpp"
// fastest_subfield_inv の twiddle 構築まわりを最適化:
//
//   1. Artin-Schreier 連鎖 c[0]=2, c[k+1]=c[k]²+c[k] (長さ 64) を **constexpr 化**。
//      PCLMUL は constexpr 不可なので 4-bit windowed clmul ベースの mul_ce / sq_ce
//      を自前実装し、 .rodata に焼き込む。
//      → init() の chain 生成段 (64 sq+xor) が完全消滅。
//
//   2. a[i] の F_2-線形結合構築を **Gray code 化**:
//        旧:  na[j] = ⊕_{k : j の k bit on} b[k]   (1 entry あたり i 回 XOR)
//        新:  na[j] = na[j ^ (1<<ctz(j))] + b[ctz(j)]   (1 entry あたり 1 回 XOR)
//      → 各 level i の構築コストが i 倍 → 1 倍に。 init の支配的部分が i 倍速化。
//
//   3. 初期 b は CHAIN の末尾 n 個を直接コピー (vector の生成 + erase を回避)。
//
// その他は fastest_subfield_inv.hpp と同一 (subfield-split inv + mul2 vectorized
// reduction + VPCLMULQDQ butterfly)。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"
#include "../../gf2-64/_shared/frob.hpp"

namespace conv_f2_64_subfield_inv_constchain {

using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// constexpr GF(2^64) mul / sq (PCLMUL は constexpr 不可なため 4-bit windowed clmul)
// 既約多項式: P(x) = x^64 + x^4 + x^3 + x + 1 (0x1B + x^64)
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
// Artin-Schreier 連鎖を compile-time 計算: c[k] = P^k(2), P(x) = x²+x。
// この多項式 / 初期値で実測すると P^63(2) = 0、 つまり c[0..62] (63 個) が
// 非零の F_2-独立な基底。 ここでは 63 個だけ持つ (元 init で pop_back していた
// 末尾 0 は最初から含めない)。 64 ではない点に注意。
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();

// =============================================================================
// subfield-split inv 関連 (fastest_subfield_inv.hpp と同一)
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

// =============================================================================
// gf2 wrapper
// =============================================================================
struct gf2 {
 u64 v;
 gf2(): v(0) {}
 gf2(u64 x): v(x) {}
 gf2 operator+(const gf2& r) const { return gf2(v ^ r.v); }
 gf2 operator-(const gf2& r) const { return gf2(v ^ r.v); }
 gf2& operator+=(const gf2& r) { v^= r.v; return *this; }
 gf2& operator-=(const gf2& r) { v^= r.v; return *this; }
 PCLMUL gf2& operator*=(const gf2& r) { v= mul(v, r.v); return *this; }
 PCLMUL gf2 operator*(const gf2& r) const { return gf2(mul(v, r.v)); }
 VPCLMUL gf2 inverse() const { return gf2(inv(v)); }
 bool operator==(const gf2& r) const { return v == r.v; }
};

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// Twiddle 表: CHAIN を起点に descent + Gray code build。
// =============================================================================
struct nim_fft_data {
 std::vector<std::vector<gf2>> a;
 VPCLMUL void init(int n) {
  if((int)a.size() > n) return;
  a.resize(n + 1);
  // 初期 basis: CHAIN の末尾 n 個 (= [c_{64-n}, ..., c_{63}])
  std::vector<gf2> b(n);
  for(int k= 0; k < n; ++k) b[k]= gf2(CHAIN[CHAIN_LEN - n + k]);
  for(int i= n;; --i) {
   if(i > 0) {
    gf2 inv_b= b.back().inverse();
    for(gf2& x : b) x*= inv_b;
   }
   // a[i] を Gray code で構築 (旧の二重 loop を排除)
   auto& na= a[i];
   na.resize(1 << i);
   na[0]= gf2(0);
   for(int j= 1; j < (1 << i); ++j) {
    int k= __builtin_ctz(j);
    na[j]= na[j ^ (1 << k)] + b[k];
   }
   if(i == 0) break;
   b.pop_back();
   for(gf2& x : b) x= x * x + x;  // P(x) = x²+x
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// nim FFT / IFFT (fastest_subfield_inv.hpp と同一: VPCLMUL + mul2 並列 reduction)
// =============================================================================
VPCLMUL inline void nim_fft(std::vector<gf2>& f) {
 int n= (int)f.size();
 std::vector<gf2> f2(n);
 int len= n;
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

}  // namespace conv_f2_64_subfield_inv_constchain

struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_subfield_inv_constchain;
  std::vector<gf2> a(n), b(m);
  for(int i= 0; i < n; ++i) a[i]= gf2(a_in[i]);
  for(int i= 0; i < m; ++i) b[i]= gf2(b_in[i]);
  auto c= nim_convolution(std::move(a), std::move(b));
  std::vector<u64> out(c.size());
  for(size_t k= 0; k < c.size(); ++k) out[k]= c[k].v;
  return out;
 }
};
