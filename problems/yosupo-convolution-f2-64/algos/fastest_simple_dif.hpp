#pragma once
#include "_common.hpp"
// fastest_simple の DIF (decimation-in-frequency) 形 baseline。
//
// 構造 (forward):
//   1. Phase A (DIT-style with shuffles): monomial → LCH 基底 (natural 順)
//   2. bit-reverse permute: LCH natural → LCH bit-reverse 順
//   3. Phase B (DIF top-down butterflies, block 単位 single ska):
//      LCH bit-reverse → eval bit-reverse 順
//
// Inverse は順序逆。
//
// LCH 基底 (Cantor): e_j(x) = ∏_{l: j_l=1} s_l(x)、 s_l(β_l)=1 (Cantor)
//   評価点 α(p) = Σ_{l: p_l=1} β_l
//
// DIF butterfly:
//   block j, level k (unit=2^k) で ska_j = s_{k-1}(α(j*2^k))
//   Cantor 基底下: s_{k-1}(β_{k+L}) = β_{L+1} (level 非依存)
//   → master M[j] = Σ_{L: j_L=1} β_{L+1} (size 2^(d-1)) で全 level の ska を提供
//   → contiguous master access が DIF の自然な benefit
//
// 注: 本 baseline では Phase A は既存 DIT のものを流用、 explicit な bit_reverse
//     を 1 回挟む。 「true DIF」 (bit-reverse を Phase A に折り込む) は次段階で検討。
//     ここでの DIF 化の主眼は butterfly 側 (DIT の per-i twiddle → DIF の
//     per-block single ska) + master 表 contiguous 化。
//
// 必要拡張: PCLMUL のみ。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_simple_dif {

using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// constexpr GF(2^64) mul / sq
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
// Cantor chain。 β_l = CHAIN[62-l]、 β_l^2 + β_l = β_{l-1}、 β_0 = 1
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
// DIF master twiddle (single contiguous table)
//   M[j] = Σ_{L: j_L=1} β_{L+1}   (j ∈ [0, 2^(d-1)))
//   level k の block j に対する ska = M[j] (level 非依存)
// =============================================================================
struct nim_fft_data {
 std::vector<u64> master;
 int max_d= -1;
 void init(int d) {
  if(d <= max_d) return;
  max_d= d;
  int sz= (d <= 1) ? 1 : (1 << (d - 1));
  master.assign(sz, 0);
  for(int j= 1; j < sz; ++j) {
   int L= __builtin_ctz(j);
   master[j]= master[j ^ (1 << L)] ^ CHAIN[61 - L];  // β_{L+1} = CHAIN[61-L]
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// in-place bit-reverse permutation
// =============================================================================
inline void bit_reverse(std::vector<u64>& f) {
 int s= (int)f.size();
 if(s <= 1) return;
 for(int i= 1, j= 0; i < s; ++i) {
  int b= s >> 1;
  while(j & b) {
   j^= b;
   b>>= 1;
  }
  j^= b;
  if(i < j) std::swap(f[i], f[j]);
 }
}

// =============================================================================
// Phase A (DIT-style with shuffles): monomial → LCH natural
//   既存 fastest_simple のものをそのまま流用。
// =============================================================================
inline void phase_a(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
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
}

// inverse Phase A (bc_to_mono、 既存 fastest_simple の逆 Phase A 流用)
inline void inv_phase_a(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= 1;
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

// =============================================================================
// Phase B: DIF top-down butterflies (single ska per block, contiguous master)
// =============================================================================
inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 phase_a(f);
 bit_reverse(f);
 const u64* M= nim_data.master.data();
 for(int k= d; k >= 1; --k) {
  int unit= 1 << k;
  int half= unit >> 1;
  int num_blocks= n / unit;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 ska= M[j];
   if(ska == 0) {
    for(int i= 0; i < half; ++i) base[half + i]^= base[i];
   } else {
    for(int i= 0; i < half; ++i) {
     u64 r= mul(base[half + i], ska);
     base[i]^= r;
     base[half + i]^= base[i];
    }
   }
  }
 }
}

inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 const u64* M= nim_data.master.data();
 for(int k= 1; k <= d; ++k) {
  int unit= 1 << k;
  int half= unit >> 1;
  int num_blocks= n / unit;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 ska= M[j];
   if(ska == 0) {
    for(int i= 0; i < half; ++i) base[half + i]^= base[i];
   } else {
    for(int i= 0; i < half; ++i) {
     base[half + i]^= base[i];
     base[i]^= mul(base[half + i], ska);
    }
   }
  }
 }
 bit_reverse(f);
 inv_phase_a(f);
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

}  // namespace conv_f2_64_simple_dif

struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_simple_dif;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
