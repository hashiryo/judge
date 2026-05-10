#pragma once
#include "_common.hpp"
// fastest_simple + **bit-reverse 中間表現** (placeholder 版):
//
// 構造的には fastest_simple とほぼ同じ。 ただし FFT 出力 / IFFT 入力に
// bit-reverse permutation を 1 回ずつ挟んで、 中間表現が bit-reverse 順に
// なるようにしてある。
//
// 注意: これは「DIF nim FFT に向けた構造の placeholder」 であって、 hot loop
// での access pattern 改善はしていない (master 表共有とも組み合わせていない)。
// bit-reverse permute 3 回ぶんのコストが乗るので simple より遅くなる。
//
// この上で、 真の DIF nim FFT (= bit-reverse intermediate を Phase A 内に
// 自然に埋め込み、 master 表 contiguous access を実現する) を実装するのが次の段階。
//
// 必要拡張: PCLMUL のみ。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_simple_bitrev {

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
// Twiddle: full table (size 2^i), 自然順
//   level i での basis = CHAIN[CHAIN_LEN-i .. CHAIN_LEN-1] (suffix 直接参照)
// =============================================================================
struct nim_fft_data {
 std::vector<std::vector<u64>> a;
 void init(int n) {
  if((int)a.size() > n) return;
  a.resize(n + 1);
  for(int i= 0; i <= n; ++i) {
   auto& na= a[i];
   na.assign(1 << i, 0);  // 全長 2^i
   const u64* basis= &CHAIN[CHAIN_LEN - i];
   for(int j= 1; j < (1 << i); ++j) {
    int k= __builtin_ctz(j);
    na[j]= na[j ^ (1 << k)] ^ basis[k];
   }
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
 // Phase B: 各 butterfly で 2 mul (lower 用 g_lo[i], upper 用 g_hi[i])
 while(len < n) {
  len*= 2;
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g_lo= g.data();
   const u64* g_hi= g.data() + half;
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    f0[i]= lo ^ mul(hi, g_lo[i]);
    f1[i]= lo ^ mul(hi, g_hi[i]);
   }
  }
 }
}

inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // 逆 Phase B
 while(len > 1) {
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g_lo= g.data();
   const u64* g_hi= g.data() + half;
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   for(int i= 0; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    f0[i]= mul(lo, g_hi[i]) ^ mul(hi, g_lo[i]);
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

// in-place bit-reverse permutation (0..s-1 で適用、 s = 2^bits)
inline void bit_reverse(std::vector<u64>& f) {
 int s= (int)f.size();
 int bits= msb(s);
 // ビット反転による in-place swap (swap pairs (i, bitrev(i)) where i < bitrev(i))
 for(int i= 1, j= 0; i < s; ++i) {
  // increment j by 1 in bit-reversed sense (Gray-code-like)
  int b= s >> 1;
  while(j & b) {
   j^= b;
   b>>= 1;
  }
  j^= b;
  if(i < j) std::swap(f[i], f[j]);
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
 // 中間表現を bit-reverse 順に切り替え (placeholder; 真の DIF ならここで permute せず
 // Phase B 内部で自然に bit-reverse 順を扱う)
 bit_reverse(f);
 bit_reverse(g);
 for(int i= 0; i < s; ++i) f[i]= mul(f[i], g[i]);
 // IFFT に渡す前に元の順序に戻す
 bit_reverse(f);
 nim_ifft(f);
 f.resize(n + m - 1);
 return f;
}

}  // namespace conv_f2_64_simple_bitrev

struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_simple_bitrev;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
