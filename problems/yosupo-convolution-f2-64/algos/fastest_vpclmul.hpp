#pragma once
#include "_common.hpp"
// fastest_v2.hpp の Phase B / IFFT-twiddle で、 各 butterfly の 2 個の PCLMUL を
// **VPCLMULQDQ** (AVX-2 拡張、 _mm256_clmulepi64_epi128) で 1 命令に圧縮。
//
// VPCLMULQDQ:
//   __m256i _mm256_clmulepi64_epi128(__m256i a, __m256i b, imm8)
//   = 2 個の独立な (64-bit a_lane * 64-bit b_lane) の polynomial 積を、
//     256-bit 内で並列に実行。
//
//   (a0, a1) と (b0, b1) を pack して 1 命令で (a0*b0, a1*b1) の 128-bit×2 を得る。
//   Ice Lake (2019+) / Zen 3 (2020+) 以降の CPU で利用可。
//
// 使い方:
//   butterfly の 2 mul (hi * g0, hi * g1) は a 側が同じ hi、 b 側が異なる g。
//   pack (hi, hi) と (g0, g1) で 1 VPCLMULQDQ で両方計算。
//
// reduction:
//   結果 256-bit から (lo0, hi0, lo1, hi1) を取り出して 2 つを並列 reduce。
//   AVX-2 の logical/shift で SIMD 化、 RED テーブル 2 lookup は scalar fallback。
//
// 効果:
//   PCLMUL throughput が約 2 倍に。 hot path の per-iter 2 PCLMUL → 1 VPCLMULQDQ。
//
// 制約: x64 + AVX-2 + VPCLMULQDQ feature 必須。 ARM / 古い x64 では fastest_v2
//       にフォールバックすべきだが、 簡単のためここでは要求のみ。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,pclmul,vpclmulqdq")
#include <immintrin.h>
#define HAS_VPCLMUL 1
#define HAS_PCLMUL 1
#elif defined(USE_SIMDE)
#define HAS_VPCLMUL 0
#define HAS_PCLMUL 1
#else
#define HAS_VPCLMUL 0
#define HAS_PCLMUL 0
#endif

namespace conv_f2_64_fastest_vpclmul {

inline const std::array<u64, 16> RED = [] {
 std::array<u64, 16> r{};
 for (int q = 0; q < 16; ++q) {
  u64 o = (u64) q ^ ((u64) q >> 1) ^ ((u64) q >> 3);
  r[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
 }
 return r;
}();

#if HAS_VPCLMUL
#define GF2_TARGET [[gnu::target("avx2,pclmul,vpclmulqdq")]]
#else
#define GF2_TARGET
#endif

struct gf2 {
 u64 v;
 gf2() : v(0) {}
 gf2(u64 x) : v(x) {}
 gf2 operator+(const gf2& r) const { return gf2(v ^ r.v); }
 gf2 operator-(const gf2& r) const { return gf2(v ^ r.v); }
 gf2& operator+=(const gf2& r) { v ^= r.v; return *this; }
 gf2& operator-=(const gf2& r) { v ^= r.v; return *this; }
#if HAS_PCLMUL
 GF2_TARGET gf2& operator*=(const gf2& r) {
  __m128i av{(long long) v, 0};
  __m128i bv{(long long) r.v, 0};
  __m128i m = _mm_clmulepi64_si128(av, bv, 0);
  u64 lo = (u64) m[0], hi = (u64) m[1];
  lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
  v = lo ^ RED[hi >> 60];
  return *this;
 }
#else
 gf2& operator*=(const gf2& r) {
  u64 a = v, b = r.v;
  u64 lo = 0, hi = 0;
  for (int i = 0; i < 64; ++i) {
   if ((b >> i) & 1) { lo ^= a << i; if (i) hi ^= a >> (64 - i); }
  }
  lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
  v = lo ^ RED[hi >> 60];
  return *this;
 }
#endif
 gf2 operator*(const gf2& r) const { gf2 t = *this; t *= r; return t; }
 gf2 pow(u64 k) const {
  gf2 res(1), a = *this;
  while (k) { if (k & 1) res *= a; a *= a; k >>= 1; }
  return res;
 }
 gf2 inv() const { return pow(~u64(1)); }
 bool operator==(const gf2& r) const { return v == r.v; }
};

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

#if HAS_VPCLMUL
// (a0*b0, a1*b1) を VPCLMULQDQ 1 命令で計算 + 並列 reduction。
GF2_TARGET inline void mul_pair(u64 a0, u64 a1, u64 b0, u64 b1, u64& out0, u64& out1) {
 // 256-bit に (a0, a1) と (b0, b1) を pack (lane 0 = a0*b0, lane 1 = a1*b1)
 __m256i av = _mm256_set_epi64x(0, a1, 0, a0);
 __m256i bv = _mm256_set_epi64x(0, b1, 0, b0);
 __m256i m = _mm256_clmulepi64_epi128(av, bv, 0);
 u64 lo0 = (u64) _mm256_extract_epi64(m, 0);
 u64 hp0 = (u64) _mm256_extract_epi64(m, 1);
 u64 lo1 = (u64) _mm256_extract_epi64(m, 2);
 u64 hp1 = (u64) _mm256_extract_epi64(m, 3);
 lo0 ^= hp0 ^ (hp0 << 1) ^ (hp0 << 3) ^ (hp0 << 4);
 lo1 ^= hp1 ^ (hp1 << 1) ^ (hp1 << 3) ^ (hp1 << 4);
 out0 = lo0 ^ RED[hp0 >> 60];
 out1 = lo1 ^ RED[hp1 >> 60];
}
#endif

struct nim_fft_data {
 std::vector<std::vector<gf2>> a;
 GF2_TARGET void init(int n) {
  if ((int) a.size() > n) return;
  a.resize(n + 1);
  std::vector<gf2> b;
  b.push_back(gf2(2));
  while (true) {
   gf2 x = b.back();
   if (x.v == 0) break;
   b.push_back(x * x + x);
  }
  b.pop_back();
  b.erase(b.begin(), b.end() - n);
  for (int i = n;; --i) {
   std::vector<gf2>& na = a[i];
   na.resize(1 << i);
   gf2 inv_b = b.back().inv();
   for (gf2& x : b) x *= inv_b;
   for (int j = 0; j < (1 << i); ++j) {
    for (int k = 0; k < i; ++k) if (j >> k & 1) na[j] += b[k];
   }
   if (i == 0) break;
   b.pop_back();
   for (gf2& x : b) x = x * x + x;
  }
 }
};
inline nim_fft_data nim_data;

GF2_TARGET inline void nim_fft(std::vector<gf2>& f) {
 int n = (int) f.size();
 std::vector<gf2> f2(n);
 int len = n;
 while (len > 1) {
  for (int l = 0; l < n; l += len) {
   for (int m = len / 4; m >= 1; m /= 2) {
    for (int t = 0; t < len; t += m * 4) {
     for (int i = 0; i < m; ++i) {
      gf2 b = f[l + t + m + i], c = f[l + t + m * 2 + i], d = f[l + t + m * 3 + i];
      f[l + t + m + i] = b + c + d;
      f[l + t + m * 2 + i] = c + d;
     }
    }
   }
   for (int i = 0; i < len / 2; ++i) {
    f2[l + i] = f[l + i * 2];
    f2[l + i + len / 2] = f[l + i * 2 + 1];
   }
  }
  std::swap(f, f2);
  len /= 2;
 }
 // Phase B: in-place + VPCLMULQDQ で 2 mul を 1 命令に圧縮
 while (len < n) {
  len *= 2;
  const std::vector<gf2>& g = nim_data.a[msb(len)];
  for (int l = 0; l < n; l += len) {
   const int half = len / 2;
   const gf2* g0 = g.data();
   const gf2* g1 = g.data() + half;
   gf2* f0 = f.data() + l;
   gf2* f1 = f.data() + l + half;
   for (int i = 0; i < half; ++i) {
    gf2 lo = f0[i], hi = f1[i];
#if HAS_VPCLMUL
    // (hi * g0, hi * g1) を 1 VPCLMULQDQ で
    u64 p0, p1;
    mul_pair(hi.v, hi.v, g0[i].v, g1[i].v, p0, p1);
    f0[i] = gf2(lo.v ^ p0);
    f1[i] = gf2(lo.v ^ p1);
#else
    f0[i] = lo + hi * g0[i];
    f1[i] = lo + hi * g1[i];
#endif
   }
  }
 }
}

GF2_TARGET inline void nim_ifft(std::vector<gf2>& f) {
 int n = (int) f.size();
 std::vector<gf2> f2(n);
 int len = n;
 while (len > 1) {
  const std::vector<gf2>& g = nim_data.a[msb(len)];
  for (int l = 0; l < n; l += len) {
   const int half = len / 2;
   const gf2* g0 = g.data();
   const gf2* g1 = g.data() + half;
   gf2* f0 = f.data() + l;
   gf2* f1 = f.data() + l + half;
   for (int i = 0; i < half; ++i) {
    gf2 lo = f0[i], hi = f1[i];
#if HAS_VPCLMUL
    // (lo * g1, hi * g0) を 1 VPCLMULQDQ で
    u64 p0, p1;
    mul_pair(lo.v, hi.v, g1[i].v, g0[i].v, p0, p1);
    f0[i] = gf2(p0 ^ p1);
    f1[i] = lo + hi;
#else
    f0[i] = lo * g1[i] + hi * g0[i];
    f1[i] = lo + hi;
#endif
   }
  }
  len /= 2;
 }
 while (len < n) {
  len *= 2;
  for (int l = 0; l < n; l += len) {
   for (int i = 0; i < len / 2; ++i) {
    f2[l + i * 2] = f[l + i];
    f2[l + i * 2 + 1] = f[l + i + len / 2];
   }
  }
  std::swap(f, f2);
  for (int l = 0; l < n; l += len) {
   for (int m = 1; m <= len / 4; m *= 2) {
    for (int t = 0; t < len; t += m * 4) {
     for (int i = 0; i < m; ++i) {
      gf2 b = f[l + t + m + i], c = f[l + t + m * 2 + i], d = f[l + t + m * 3 + i];
      f[l + t + m + i] = b + c;
      f[l + t + m * 2 + i] = c + d;
     }
    }
   }
  }
 }
}

GF2_TARGET inline std::vector<gf2> nim_convolution(std::vector<gf2> f, std::vector<gf2> g) {
 int n = (int) f.size(), m = (int) g.size();
 int s = (int) ceil_pow2(u32(n + m - 1));
 f.resize(s); g.resize(s);
 nim_data.init(msb(s));
 nim_fft(f); nim_fft(g);
 for (int i = 0; i < s; ++i) f[i] *= g[i];
 nim_ifft(f);
 f.resize(n + m - 1);
 return f;
}

}  // namespace

struct Solver {
 GF2_TARGET static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_fastest_vpclmul;
  std::vector<gf2> a(n), b(m);
  for (int i = 0; i < n; ++i) a[i] = gf2(a_in[i]);
  for (int i = 0; i < m; ++i) b[i] = gf2(b_in[i]);
  auto c = nim_convolution(std::move(a), std::move(b));
  std::vector<u64> out(c.size());
  for (size_t k = 0; k < c.size(); ++k) out[k] = c[k].v;
  return out;
 }
};
