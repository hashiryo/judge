#pragma once
#include "_common.hpp"
// fastest_simple の DIF (decimation-in-frequency) 形 — true DIF 版。
// fastest_simple_dif.hpp (= placeholder、 explicit bit-reverse 版) と並列比較用。
//
// 構造 (forward):
//   1. Phase A (true DIF, no shuffle, no bit-reverse):
//      monomial 基底 → LCH 自然順 を直接生成。
//      再帰: P = Q_0 + s_{d-1}(x) · Q_1 と分解 (s_{d-1} = Cantor subspace poly)、
//      Q_0/Q_1 を半分メモリに格納し、 各半分を再帰 bc。
//      Cantor 基底下、 s_k(x) = Σ_{l ⊆ k bits} x^{2^l} (Lucas/Pascal mod 2)。
//   2. Phase B (DIF top-down butterflies, block 単位 single ska):
//      LCH natural → eval natural。
//
// Inverse は bc_to_mono = forward 操作の reverse 順 (各 XOR は self-inverse)。
//
// LCH 基底 (Cantor): e_j(x) = ∏_{l: j_l=1} s_l(x)、 s_l(β_l) = 1。
//   評価点 α(p) = Σ_{l: p_l=1} β_l。
//
// DIF butterfly twiddle:
//   block j, level k (unit=2^k) で ska_j = s_{k-1}(α(j·2^k))。
//   Cantor 帰納で s_{k-1}(β_{k+L}) = β_{L+1} (level 非依存)。
//   → master M[j] = Σ_{L: j_L=1} β_{L+1} (size 2^(d-1)) で全 level の ska を提供。
//   → contiguous master access が DIF の自然な benefit。
//
// 必要拡張: PCLMUL のみ。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_simple_dif_true {

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
// Submask 表 (proper submasks of k for k=1..MAX_SUB_IDX)
//   subs[k][.] = (1 << l) for proper submask l of k、 sub_n[k] = 個数。
//   n は int サイズ ≤ 2^31-1 なので d ≤ 30、 sub_idx ≤ 29 で MAX_SUB_IDX=30。
//   k ≤ 30 で popcount(k) max = 4 (popcount(15)=popcount(23)=popcount(30)=4 等)、
//   proper submask 数 = 2^4 - 1 = 15 → MAX_SUBS_PER_K = 16 で足りる。
// =============================================================================
constexpr int MAX_SUB_IDX= 30;
constexpr int MAX_SUBS_PER_K= 16;
struct SubmaskTable {
 array<array<int, MAX_SUBS_PER_K>, MAX_SUB_IDX + 1> subs;
 array<int, MAX_SUB_IDX + 1> sub_n;
};
constexpr SubmaskTable build_submask_tbl() {
 SubmaskTable t{};
 for(int k= 1; k <= MAX_SUB_IDX; ++k) {
  int n= 0;
  int s= k;
  while(true) {
   s= (s - 1) & k;
   t.subs[k][n++]= 1 << s;
   if(s == 0) break;
  }
  t.sub_n[k]= n;
 }
 return t;
}
constexpr auto SUBMASK_TBL= build_submask_tbl();

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
  int sub_idx= level - 1;  // s_{sub_idx} = subspace poly of dim 2^sub_idx
  const int sub_n= SUBMASK_TBL.sub_n[sub_idx];
  const int* subs= SUBMASK_TBL.subs[sub_idx].data();
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= half - 1; i >= 0; --i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    for(int idx= 0; idx < sub_n; ++idx) p[i + subs[idx]]^= q;
   }
  }
 }
}

// inverse bc (bc_to_mono): forward の XOR 列を逆順で適用。
//   各 XOR は self-inverse、 但し src/dst dependency があるので順序が重要。
//   outer level: 小→大 (forward の逆)、 block 内 i: 0 → half-1 (forward の逆)、
//   submask l: 順序は forward の逆 (proper submask を逆順で)。
inline void bc_to_mono(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= 2; level <= d; ++level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  const int sub_n= SUBMASK_TBL.sub_n[sub_idx];
  const int* subs= SUBMASK_TBL.subs[sub_idx].data();
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= 0; i < half; ++i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    // forward の逆順で適用 (subs を逆順)
    for(int idx= sub_n - 1; idx >= 0; --idx) p[i + subs[idx]]^= q;
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
 bc_to_lch(f.data(), n);
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

}  // namespace conv_f2_64_simple_dif_true

struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_simple_dif_true;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
