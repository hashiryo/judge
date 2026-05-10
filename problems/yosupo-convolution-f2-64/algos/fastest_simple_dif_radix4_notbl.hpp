#pragma once
#include "_common.hpp"
// fastest_simple_dif_radix4.hpp の **submask table 不要版**。
// (radix-4 + inline submask while ループ、 SUBMASK_TBL なし)
// 元コメント:
//
// 通常 (radix-2) DIF は level k=d, d-1, ..., 1 を 1 段ずつ降りながら butterfly する。
// 各 level で n/2 mul、 d level で 計 nd/2 mul、 メモリも 1 element / level 触る。
//
// radix-4 では level pair (k, k-1) を 1 パスにまとめる:
//   block size 2^k を 4 quarters (Q_0, Q_1, Q_2, Q_3) に分け、 4-way butterfly で
//   4 出力 (R_C00, R_C01, R_C10, R_C11) を一括生成。
//
// LCH 構造から:
//   P = Q_0 + s_{k-2}·Q_1 + s_{k-1}·Q_2 + s_{k-2}s_{k-1}·Q_3
//   評価点 γ で P(γ) = Q_0(γ) + u·Q_1(γ) + v·Q_2(γ) + uv·Q_3(γ)
//     u = s_{k-2}(γ_block), v = s_{k-1}(γ_block)
//   4 つの coset で u, v は (u, v), (u+1, v), (u+β_1, v+1), (u+β_1+1, v+1)
//   (s_{k-2}(β_{k-2}) = 1, s_{k-2}(β_{k-1}) = β_1, s_{k-1}(β_{k-1}) = 1)
//
// 2 段の radix-2 として因数分解:
//   Step 1 (twiddle v): A_0 = Q_0 ^ v Q_2; A_2 = A_0 ^ Q_2;
//                       A_1 = Q_1 ^ v Q_3; A_3 = A_1 ^ Q_3
//   Step 2 lower (twiddle u_lo = u): R_C00 = A_0 ^ u A_1; R_C01 = R_C00 ^ A_1
//   Step 2 upper (twiddle u_hi = u + β_1): R_C10 = A_2 ^ u_hi A_3; R_C11 = R_C10 ^ A_3
//
// mul 数: 4 per group of 4 elements = radix-2 と同じ (per-element 1 mul/double-level)。
// 主な利点はメモリトラフィック半減 (radix-2 では 2 level で 各要素を 2 回 read/write、
// radix-4 では 1 回ずつ)。
//
// twiddle: master M[j] = Σ_{L: j_L=1} β_{L+1} は level 非依存なので
//   v = M[j_k]、 u_lo = M[2 j_k]、 u_hi = M[2 j_k + 1] で全部賄える。
//
// 必要拡張: PCLMUL のみ。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_simple_dif_radix4_notbl {

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
// Phase A (true DIF, no shuffle): monomial → LCH natural
//   level len = n, n/2, ..., 4 (top-down)。 各 level で contiguous block ごとに
//   "polynomial division by s_{k-1}" (k = log2(len)) の feedback XOR を流す。
//   その後 block を半分に分けて次 level (smaller len)。
//
//   Cantor で s_k(x) = Σ_{l: l⊆k bits} x^{2^l}。
//   block 内 division は top-down (i = half-1 → 0):
//     for each l ∈ S_{k-1} \ {k-1} (= proper submask of k-1):
//       poly[base + i + 2^l] ^= poly[base + i + half]
//   Q_1 (高位 half) はそのまま位置 [half, len) に残る。
// =============================================================================
inline void bc_to_lch(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= d; level >= 2; --level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= half - 1; i >= 0; --i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    int s= sub_idx;
    while(true) {
     s= (s - 1) & sub_idx;
     p[i + (1 << s)]^= q;
     if(s == 0) break;
    }
   }
  }
 }
}

inline void bc_to_mono(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= 2; level <= d; ++level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= 0; i < half; ++i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    // submask 順序は inner で commutative
    int s= sub_idx;
    while(true) {
     s= (s - 1) & sub_idx;
     p[i + (1 << s)]^= q;
     if(s == 0) break;
    }
   }
  }
 }
}

// =============================================================================
// Phase B: DIF top-down radix-4 butterflies
//   level pair (k, k-1) を 1 パスで処理。 d 奇数なら最後に level 1 の radix-2 1 段。
// =============================================================================
inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 bc_to_lch(f.data(), n);
 const u64* M= nim_data.master.data();
 int k= d;
 while(k >= 2) {
  int unit= 1 << k;
  int q= unit >> 2;       // 2^{k-2}
  int q2= q << 1;          // 2^{k-1}
  int q3= q + q2;          // 3·2^{k-2}
  int num_blocks= n / unit;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 v= M[j];
   u64 u_lo= M[2 * j];
   u64 u_hi= M[2 * j + 1];  // u_hi = u_lo ^ β_1
   for(int i= 0; i < q; ++i) {
    u64 Q0= base[i];
    u64 Q1= base[q + i];
    u64 Q2= base[q2 + i];
    u64 Q3= base[q3 + i];
    // Step 1: twiddle v
    u64 A0, A1, A2, A3;
    if(v == 0) {
     A0= Q0;
     A1= Q1;
    } else {
     A0= Q0 ^ mul(Q2, v);
     A1= Q1 ^ mul(Q3, v);
    }
    A2= A0 ^ Q2;
    A3= A1 ^ Q3;
    // Step 2 lower: twiddle u_lo
    u64 R0, R1;
    if(u_lo == 0) {
     R0= A0;
    } else {
     R0= A0 ^ mul(A1, u_lo);
    }
    R1= R0 ^ A1;
    // Step 2 upper: twiddle u_hi (常に non-zero、 u_lo ^ β_1)
    u64 R2= A2 ^ mul(A3, u_hi);
    u64 R3= R2 ^ A3;
    base[i]= R0;
    base[q + i]= R1;
    base[q2 + i]= R2;
    base[q3 + i]= R3;
   }
  }
  k-= 2;
 }
 // d 奇数なら最後に level 1 の radix-2 (block size 2)
 if(k == 1) {
  int num_blocks= n >> 1;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * 2;
   u64 ska= M[j];
   if(ska == 0) {
    base[1]^= base[0];
   } else {
    u64 r= mul(base[1], ska);
    base[0]^= r;
    base[1]^= base[0];
   }
  }
 }
}

inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 const u64* M= nim_data.master.data();
 // 逆順: forward が pair (d, d-1), (d-2, d-3), ..., の順で top-down だったので
 // inverse は逆 (smallest pair_level → d)。 d 奇数なら先頭で level 1 radix-2 を invert。
 int k_start;
 if(d & 1) {
  // level 1 radix-2 inverse 先行
  int num_blocks= n >> 1;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * 2;
   u64 ska= M[j];
   if(ska == 0) {
    base[1]^= base[0];
   } else {
    base[1]^= base[0];
    base[0]^= mul(base[1], ska);
   }
  }
  k_start= 3;  // 次の pair の upper level
 } else {
  k_start= 2;
 }
 for(int level_hi= k_start; level_hi <= d; level_hi+= 2) {
  // pair (level_hi, level_hi - 1) inverse、 block size = 2^level_hi
  int unit= 1 << level_hi;
  int q= unit >> 2;
  int q2= q << 1;
  int q3= q + q2;
  int num_blocks= n / unit;
  for(int j= 0; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 v= M[j];
   u64 u_lo= M[2 * j];
   u64 u_hi= M[2 * j + 1];
   for(int i= 0; i < q; ++i) {
    u64 R0= base[i];
    u64 R1= base[q + i];
    u64 R2= base[q2 + i];
    u64 R3= base[q3 + i];
    // Step 2 inverse upper: A_3 = R2 ^ R3, A_2 = R2 ^ u_hi·A_3
    u64 A3= R2 ^ R3;
    u64 A2= R2 ^ mul(A3, u_hi);
    // Step 2 inverse lower: A_1 = R0 ^ R1, A_0 = R0 ^ u_lo·A_1
    u64 A1= R0 ^ R1;
    u64 A0;
    if(u_lo == 0) A0= R0;
    else A0= R0 ^ mul(A1, u_lo);
    // Step 1 inverse: Q_2 = A_0 ^ A_2, Q_0 = A_0 ^ v·Q_2; Q_3 = A_1 ^ A_3, Q_1 = A_1 ^ v·Q_3
    u64 Q2= A0 ^ A2;
    u64 Q3= A1 ^ A3;
    u64 Q0, Q1;
    if(v == 0) {
     Q0= A0;
     Q1= A1;
    } else {
     Q0= A0 ^ mul(Q2, v);
     Q1= A1 ^ mul(Q3, v);
    }
    base[i]= Q0;
    base[q + i]= Q1;
    base[q2 + i]= Q2;
    base[q3 + i]= Q3;
   }
  }
 }
 bc_to_mono(f.data(), n);
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

}  // namespace conv_f2_64_simple_dif_radix4_notbl

struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_simple_dif_radix4_notbl;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
