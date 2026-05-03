#pragma once
#include "_common.hpp"
// fastest.hpp の Phase B (twiddle mul) を:
//   1. in-place 化 (f2 への copy + swap を排除、 pair 単位で lo, hi を register に
//      load してから上書き)
//   2. 2× unroll (i, i+1 の butterfly を交互発行、 PCLMUL を 4 個 in-flight)
// で改善。 fft の Phase B と ifft 全部に適用。 Phase A の deinterleave は f2 必要
// なので元のまま。
//
// 期待: PCLMUL latency ~5 cyc を OoO で隠蔽、 hot path で ~30-50% 高速化見込み
//
// アルゴリズム: 標数 2 の有限体 F_{2^64} 上の additive FFT (a.k.a. nim FFT)
//   - 基底列 b: b[0] = 2 ∈ F_{2^64}、 b[i+1] = b[i]^2 + b[i] (Artin-Schreier)
//     これらは独立で、 b[0..63] が F_{2^64} の F_2 上の基底になる
//   - n = 2^k 点の FFT は: 多項式 f(X) を線形 polylog で部分評価
//   - precompute: 各 level i ∈ [0, k) について 2^i 個の twiddle factor を作る
//     (b[0..i-1] の F_2 線形結合)
//   - fft: 2 段階
//     Phase A (butterfly): f[t+m+i], f[t+2m+i], f[t+3m+i] → 加算で前処理
//     Phase B (twiddle mul): 各 level で twiddle と XOR で点評価を伸長
//   - ifft: 逆順
//   - convolution: pad → fft(a), fft(b) → 点ごと積 → ifft → trim
//
// 計算量: O((n+m) log(n+m)) F_{2^64} 演算
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,sse4.2")
#include <immintrin.h>
#endif

namespace conv_f2_64_fastest_v2 {

inline const std::array<u64, 16> RED = [] {
 std::array<u64, 16> r{};
 for (int q = 0; q < 16; ++q) {
  u64 o = (u64) q ^ ((u64) q >> 1) ^ ((u64) q >> 3);
  r[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
 }
 return r;
}();

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define GF2_TARGET [[gnu::target("pclmul,sse4.2")]]
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
 GF2_TARGET gf2& operator*=(const gf2& r) {
  __m128i av{(long long) v, 0};
  __m128i bv{(long long) r.v, 0};
  __m128i m = _mm_clmulepi64_si128(av, bv, 0);
  u64 lo = (u64) m[0], hi = (u64) m[1];
  lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
  v = lo ^ RED[hi >> 60];
  return *this;
 }
 GF2_TARGET gf2 operator*(const gf2& r) const { gf2 t = *this; t *= r; return t; }
 gf2 pow(u64 k) const {
  gf2 res(1), a = *this;
  while (k) { if (k & 1) res *= a; a *= a; k >>= 1; }
  return res;
 }
 gf2 inv() const { return pow(~u64(1)); }  // = pow(2^64 - 2)
 bool operator==(const gf2& r) const { return v == r.v; }
};

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// Twiddle 表: a[i] は level i (= 2^i 点 FFT) 用の 2^i サイズ vector
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
  // b の末尾 n 個だけ残す
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
 // Phase A: 連続的な butterfly + bit-reverse 風置換、 length が小さくなる方向
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
 // Phase B: 各 level で twiddle 乗算 (in-place + 2× unroll)
 while (len < n) {
  len *= 2;
  const std::vector<gf2>& g = nim_data.a[msb(len)];
  for (int l = 0; l < n; l += len) {
   const int half = len / 2;
   const gf2* g0 = g.data();
   const gf2* g1 = g.data() + half;
   gf2* f0 = f.data() + l;
   gf2* f1 = f.data() + l + half;
   int i = 0;
   for (; i + 1 < half; i += 2) {
    gf2 lo0 = f0[i],   hi0 = f1[i];
    gf2 lo1 = f0[i+1], hi1 = f1[i+1];
    gf2 hi0_g0 = hi0 * g0[i];
    gf2 hi0_g1 = hi0 * g1[i];
    gf2 hi1_g0 = hi1 * g0[i+1];
    gf2 hi1_g1 = hi1 * g1[i+1];
    f0[i]   = lo0 + hi0_g0;
    f1[i]   = lo0 + hi0_g1;
    f0[i+1] = lo1 + hi1_g0;
    f1[i+1] = lo1 + hi1_g1;
   }
   for (; i < half; ++i) {
    gf2 lo = f0[i], hi = f1[i];
    f0[i] = lo + hi * g0[i];
    f1[i] = lo + hi * g1[i];
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
   int i = 0;
   for (; i + 1 < half; i += 2) {
    gf2 lo0 = f0[i],   hi0 = f1[i];
    gf2 lo1 = f0[i+1], hi1 = f1[i+1];
    gf2 lo0_g1 = lo0 * g1[i];
    gf2 hi0_g0 = hi0 * g0[i];
    gf2 lo1_g1 = lo1 * g1[i+1];
    gf2 hi1_g0 = hi1 * g0[i+1];
    f0[i]   = lo0_g1 + hi0_g0;
    f1[i]   = lo0 + hi0;
    f0[i+1] = lo1_g1 + hi1_g0;
    f1[i+1] = lo1 + hi1;
   }
   for (; i < half; ++i) {
    gf2 lo = f0[i], hi = f1[i];
    f0[i] = lo * g1[i] + hi * g0[i];
    f1[i] = lo + hi;
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
 GF2_TARGET static std::string run(const std::string& input) {
  using namespace conv_f2_64_fastest_v2;
  std::istringstream in(input);
  std::ostringstream out;
  int n, m;
  in >> n >> m;
  std::vector<gf2> a(n), b(m);
  for (auto& x : a) { u64 v; in >> v; x = gf2(v); }
  for (auto& x : b) { u64 v; in >> v; x = gf2(v); }
  auto c = nim_convolution(std::move(a), std::move(b));
  for (int k = 0; k < (int) c.size(); ++k) {
   out << c[k].v;
   out << (k + 1 == (int) c.size() ? '\n' : ' ');
  }
  return std::move(out).str();
 }
};
