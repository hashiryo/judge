#pragma once
#include "_common.hpp"
// simple + **master twiddle 表共有** (Step 2-α):
//
// 観察: a[i][j] = a[n][j << (n-i)] (n = max level)。 つまり level i の表は
// master = a[n] の stride 2^(n-i) sub-sample。 全 level の twiddle 情報は master
// だけで足りる。
//
// 実装:
//   - init() は master (size 2^n) のみ構築。 per-level の dense table は持たない。
//   - hot loop は master[i << (n-level)] を strided access。
//
// 効果:
//   - twiddle メモリが 2^(n+1) - 1 → 2^n entry (約半減)。 s=2^20 で 16MB → 8MB。
//   - 副作用: hot loop のアクセスが contiguous → strided になり、 prefetcher 依存。
//     特に深い level (k 小) で stride が大きく、 cache miss 増の可能性。
//
// なお lower/upper の関係: 任意 level k で
//   g_lo[i] = master[i * stride]
//   g_hi[i] = master[i * stride + 2^(n-1)]   (level に無関係に一定オフセット 2^(n-1))
// なので g_lo / g_hi の base ポインタは固定で良い。
//
// 必要拡張: PCLMUL のみ。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_simple_master {

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
// c[62] = 1, c[63] = 0 なので非零 basis は c[0..62] の 63 個。
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// Twiddle: master 表のみ (size 2^max_level)、 全 level がここから stride で参照
// =============================================================================
struct nim_fft_data {
 std::vector<u64> master;
 int max_level= -1;
 void init(int n) {
  if(n <= max_level) return;
  max_level= n;
  master.assign(1ULL << n, 0);
  // Gray code build (level n の full table)
  const u64* basis= &CHAIN[CHAIN_LEN - n];
  for(int j= 1; j < (1 << n); ++j) {
   int k= __builtin_ctz(j);
   master[j]= master[j ^ (1 << k)] ^ basis[k];
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// nim FFT (out-of-place、 scalar mul、 全 twiddle 直接参照)
// =============================================================================
inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // Phase A: pure XOR butterfly + shuffle (deinterleave)
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
 // Phase B: master を stride access。 g_lo[i]=master[i*stride], g_hi[i]=master[i*stride + 2^(N-1)]
 const int N= nim_data.max_level;
 const u64* m_lo= nim_data.master.data();
 const u64* m_hi= nim_data.master.data() + (1LL << (N - 1));
 while(len < n) {
  len*= 2;
  const int level= msb(len);
  const int stride= 1 << (N - level);
  const int half= len / 2;
  for(int l= 0; l < n; l+= len) {
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    int g_idx= i * stride;
    f0[i]= lo ^ mul(hi, m_lo[g_idx]);
    f1[i]= lo ^ mul(hi, m_hi[g_idx]);
   }
  }
 }
}

inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // 逆 Phase B
 const int N= nim_data.max_level;
 const u64* m_lo= nim_data.master.data();
 const u64* m_hi= nim_data.master.data() + (1LL << (N - 1));
 while(len > 1) {
  const int level= msb(len);
  const int stride= 1 << (N - level);
  const int half= len / 2;
  for(int l= 0; l < n; l+= len) {
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    int g_idx= i * stride;
    f0[i]= mul(lo, m_hi[g_idx]) ^ mul(hi, m_lo[g_idx]);
    f1[i]= lo ^ hi;
   }
  }
  len/= 2;
 }
 // 逆 Phase A: deinterleave + 逆 butterfly
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

inline std::vector<u64> nim_convolution(std::vector<u64> f, std::vector<u64> g) {
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

}  // namespace conv_f2_64_simple_master

struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_simple_master;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
